import re
import os

from json import JSONEncoder
from .lc_material import LCMaterial
import bpm_backend as bpm

import dtmm
dtmm.conf.set_verbose(2)

from vtk import vtkImageData, vtkXMLImageDataReader, vtkXMLImageDataWriter
from vtk.util import numpy_support as vn

import numpy as np

import multiprocessing
import pyfftw 

from warnings import simplefilter
simplefilter(action='ignore', category=FutureWarning)


class LightPropagator:
    """The LightPropagator class allows to propagate optical fields through a LC sample.
    The simulation and choice of backend is done by calling the method ``propagate_field``. 
    For each wavelength specified at construction, two simulations are done: one with a light
    source polarised along x, and one with a light source polarised along y. This allows us to
    fully caracterize the transmission of the LC sample and reconstruct any kind of optical
    micrograph.

    Currently, we only support single plane-wave source with input fields propagating along
    ``z``. If you want to take into account input numerical aperture effect, you need to use
    the ``dtmm`` backend directly (this may change in a future version). However, we do
    support numerical aperture effect due to the microscope objective (i.e. on the output
    side).

    Parameters
    ----------
    material : :class:`~nemaktis.lc_material.LCMaterial` object
    wavelengths : array-like object
        An array containing all the wavelengths of the spectrum for the light source.
    numerical_aperture : float
        Sets the numerical aperture for the microscope objective.
    """
    def __init__(self, *, material, wavelengths, numerical_aperture):
        if not isinstance(material, LCMaterial):
            raise TypeError("material should be a LCMaterial object")
        self._material = material
        self._numerical_aperture = numerical_aperture
        self._wavelengths = list(wavelengths)

    @property
    def material(self):
        """Returns the current LC material"""
        return self._material

    def propagate_fields(self, *, method, bulk_filename=None):
        """Propagate optical fields through the LC sample using the specified backend.

        Parameters
        ----------
        method : "bpm" | "dtmm"
            If equal to "bpm", the beam propagation backend will be used. Should be used
            if accuracy is privileged over speed.

            If equal to "dtmm", the diffractive transfer matrix backend
            will be used (in its simplest version). Should be used if
            speed is privileged over accuracy.
        bulk_filename : None or string
            If none, the backend will not export the bulk value of the optical fields in the
            LC layer.  Else, the bulk fields values will be exported to a vti file whose
            basename is set by this parameter.
        """
        if method=="bpm":
            return self._bpm_propagation(bulk_filename)
        elif method=="dtmm":
            return self._dtmm_propagation(bulk_filename)
        else:
            match = re.compile("dtmm\((\d+)\)").match(method)
            if match:
                return self._dtmm_propagation(
                    bulk_filename, diffraction=int(match.group(1)))
            else:
                raise Exception("Unrecognised method, should be 'bpm' or 'dtmm'")

    def _bpm_propagation(self, bulk_filename):
        print("{ Running beam propagation backend }\n")

        director_field = self._material.director_field
        dims = director_field.get_mesh_dimensions()
        spacings = director_field.get_mesh_spacings()
        json_str = JSONEncoder().encode({
            "Algorithm settings": {
                "General": {
                    "LC field type":                      "Director",
                    "Results folder name":                "" },
                "Beam propagation": {
                    "N Woodbury steps":                   2,
                    "Boundary condition types":	          ["Periodic", "Periodic"],
                    "Wide angle corrections":             [False, False] }},
            "Physics settings": {
                "Initial conditions": {
                    "Beam profile":                       "UniformBeam",
                    "LC field file":                      "",
                    "Mesh dimensions":                    dims,
                    "Mesh spacings":                      spacings,
                    "Basis convention":                   "XYZ" },
                "Coefficients": {
                    "no":                                 str(self._material.no),
                    "ne":                                 str(self._material.ne),
                    "nhost":                              str(self._material.nhost),
                    "nin":                                str(self._material.nin),
                    "Wavelengths":                        self._wavelengths }},
            "Postprocessor settings": {
                "Bulk output": {
                    "Activate":	                          bulk_filename is not None,
                    "Base name":                          bulk_filename if bulk_filename is not None else ""},
                "Screen output": {
                    "Activate":	                          True,
                    "Base name":                          "",
                    "Isotropic layer thicknesses":        self._material.iso_layer_thicknesses,
                    "Isotropic layer refractive indices": self._material.iso_layer_indices,
                    "Focalisation z-shift":               0,
                    "Numerical aperture":                 self._numerical_aperture }}})

        if director_field.mask_type is None:
            director_vals = director_field.vals.flatten()
        else:
            mask_vals = np.expand_dims(2*director_field.mask_vals.astype(float)-1, 3)
            director_vals = np.concatenate((director_field.vals, mask_vals), 3).flatten()

        N_fields_vals = len(self._wavelengths)*4*dims[0]*dims[1]
        data_out = bpm.run_backend(json_str, director_vals, N_fields_vals)
        print("")

        output_fields = OpticalFields(
            wavelengths = self._wavelengths,
            mesh_lengths = (spacings[0]*(dims[0]-1), spacings[1]*(dims[1]-1)),
            mesh_dimensions = (dims[0], dims[1]))
        output_fields.vals = data_out.reshape(
            (len(self._wavelengths), 4, dims[1], dims[0]))/np.sqrt(2)
        return output_fields


    def _dtmm_propagation(self, bulk_filename, diffraction=1):
        print("{ Running diffraction transfer matrix backend }\n")

        director_field = self._material.director_field
        dims = director_field.get_mesh_dimensions()
        spacings = director_field.get_mesh_spacings()

        if np.abs(spacings[0]-spacings[1])>1e-6:
            # 2D simulation with an artificial spacings along the normal
            if spacings[0]==0:
                spacings[0] = spacings[1] 
            elif spacings[1]==0:
                spacings[1] = spacings[0] 
            else:
                raise Exception("dtmm supports only uniform spacings in the XY plane.")

        if isinstance(self._material.ne, str):
            if "lambda" in self._material.ne:
                print("Warning: dtmm does not support dispersive index; Using ne(0.6µm) instead")
            l = 0.6
            ne = eval(self._material.ne.replace("lambda","l").replace("^","**"))
        else:
            ne = self._material.ne
        if isinstance(self._material.no, str):
            if "lambda" in self._material.no:
                print("Warning: dtmm does not support dispersive index; Using no(0.6µm) instead")
            l = 0.6
            no = eval(self._material.no.replace("lambda","l").replace("^","**"))
        else:
            no = self._material.no
        if isinstance(self._material.nhost, str):
            if "lambda" in self._material.nhost:
                print("Warning: dtmm does not support dispersive index; Using nhost(0.6µm) instead")
            l = 0.6
            nhost = eval(self._material.nhost.replace("lambda","l").replace("^","**"))
        else:
            nhost = self._material.nhost
        if isinstance(self._material.nhost, str):
            if "lambda" in self._material.nin:
                print("Warning: dtmm does not support dispersive index; Using nin(0.6µm) instead")
            l = 0.6
            nin = eval(self._material.nin.replace("lambda","l").replace("^","**"))
        else:
            nin = self._material.nin

        if len(self._material.iso_layer_indices)!=0:
            print("Warning: specified isotropic layers will be ignored since this feature is "+
                  "not yet supported in dtmm.")
        print("")
        
        optical_data = dtmm.director2data(
           director_field.vals, mask = director_field.mask_vals,
           no = no, ne = ne, nhost = nhost,
           thickness = spacings[2]/spacings[0]*np.ones(dims[2]))

        wavelengths = 1000*np.array(self._wavelengths)
        field_data_in = dtmm.illumination_data(
            (dims[1],dims[0]), wavelengths, pixelsize = 1000*spacings[0], n = nin)
        field_data_out = dtmm.transfer_field(
            field_data_in, optical_data, nin=nin,
            betamax=self._numerical_aperture, diffraction=diffraction,
            ret_bulk=bulk_filename is not None)[0]
        print("")

        if bulk_filename is not None:
            print("{ Saving optical fields to "+bulk_filename+".vti }")
            lengths = director_field.get_mesh_lengths()

            vti_data = vtkImageData()
            vti_data.SetDimensions(dims[0], dims[1], dims[2]+1)
            vti_data.SetOrigin(-lengths[0]/2, -lengths[1]/2, -lengths[2]/2)
            vti_data.SetSpacing(spacings[0], spacings[1], spacings[2]*(dims[2]-1)/dims[2])

            for i in range(0,len(wavelengths)):
                E_inputX = field_data_out[:-1,0,i,[0,2],:,:].transpose(
                    (1,0,2,3)).reshape((2,dims[0]*dims[1]*(dims[2]+1))).transpose()
                E_inputY = field_data_out[:-1,1,i,[0,2],:,:].transpose(
                    (1,0,2,3)).reshape((2,dims[0]*dims[1]*(dims[2]+1))).transpose()

                E_real_inputX = vn.numpy_to_vtk(np.real(E_inputX))
                E_real_inputX.SetName("E_real_inputX_%sum" % wavelengths[i])
                vti_data.GetPointData().AddArray(E_real_inputX)
                E_imag_inputX = vn.numpy_to_vtk(np.imag(E_inputX))
                E_imag_inputX.SetName("E_imag_inputX_%sum" % wavelengths[i])
                vti_data.GetPointData().AddArray(E_imag_inputX)

                E_real_inputY = vn.numpy_to_vtk(np.real(E_inputY))
                E_real_inputY.SetName("E_real_inputY_%sum" % wavelengths[i])
                vti_data.GetPointData().AddArray(E_real_inputY)
                E_imag_inputY = vn.numpy_to_vtk(np.imag(E_inputY))
                E_imag_inputY.SetName("E_imag_inputY_%sum" % wavelengths[i])
                vti_data.GetPointData().AddArray(E_imag_inputY)

            writer = vtkXMLImageDataWriter()
            writer.SetFileName(bulk_filename+".vti")
            writer.SetInputData(vti_data)
            writer.Write()

            # We only keep the last slice to compute micrographs
            field_data_out = field_data_out[-1,:,:,:,:,:]

        output_fields = OpticalFields(
            wavelengths = self._wavelengths,
            mesh_lengths = (spacings[0]*(dims[0]-1), spacings[1]*(dims[1]-1)),
            mesh_dimensions = (dims[0], dims[1]))
        output_fields.vals = field_data_out[:,:,[0,2],:,:].transpose(
            (1,0,2,3,4)).reshape((len(wavelengths),4,dims[1],dims[0]))
        return output_fields


class OpticalFields:
    """The OpticalFields object stores the mesh information of the transverse mesh (plane mesh
    orthogonal to the z-direction, default altitude of 0) and the optical fields values on
    this mesh.  Since this python package is mainly used to reconstruct micrographs, we only
    store internally the complex horizontal electric field for two simulation: one with a
    light source polarised along ``x``, and the other with a light source polarised along
    ``y``.  In case multiple wavelengths were used in the simulation, we store these
    quantities separately for each wavelength.

    This class is initialised given either a wavelength array and the lengths and dimensions
    of the 2D mesh for the transverse fields or a path to a vti file containing previously
    calculated optical fields and mesh details.
    
    In the first version of this constructor:

    .. code-block:: python
    
        optical_fields = OpticalFields(
            wavelengths=[l0,l1,...,lN], mesh_lengths=(Lx,Ly), mesh_dimensions=(Nx,Ny))

    the actual values of the transverse fields needs to be provided later using the raw
    setter method fields_vals (shape (N_wavelengths,4,Ny,Nx)).

    In the second version of this constructor:

    .. code-block:: python
    
        optical_fields = OpticalFields(vti_file="path to vti file")

    the values of the wavelengths and transverse fields are automatically assigned from the
    vti file.
    """
    def __init__(self, **kwargs):
        if len(kwargs)==1:
            if not os.path.isfile(kwargs["vti_file"]):
                raise Exception("VTI file does not exists")

            print("{ Initializing optical fields from "+kwargs["vti_file"]+" }")

            reader = vtkXMLImageDataReader()
            reader.SetFileName(kwargs["vti_file"])
            reader.Update()

            point_data = reader.GetOutput().GetPointData()
            dims = np.array(reader.GetOutput().GetDimensions())
            spacings = np.array(reader.GetOutput().GetSpacing())

            if dims[2]!=1:
                raise Exception("The specified vti file should include 2D data")

            self._wavelengths = []
            for i in range(0,point_data.GetNumberOfArrays()):
                array_name = point_data.GetArrayName(i)
                match = re.compile("E_real_inputX_(\d+\.?\d*)um").match(array_name)
                if match:
                    self._wavelengths.append(match.group(1))
                    if not point_data.HasArray("E_imag_inputX_%sum" % match.group(1)):
                        raise Exception("Missing vtk array in the given vti file")
                    if not point_data.HasArray("E_real_inputY_%sum" % match.group(1)):
                        raise Exception("Missing vtk array in the given vti file")
                    if not point_data.HasArray("E_imag_inputY_%sum" % match.group(1)):
                        raise Exception("Missing vtk array in the given vti file")

            self._vals = pyfftw.empty_aligned(
                (len(self._wavelengths),4,dims[1],dims[0]), dtype="complex128")
            self._fft_vals = pyfftw.empty_aligned(
                (len(self._wavelengths),4,dims[1],dims[0]), dtype="complex128")
            self._fft_plan = pyfftw.FFTW(
                self._vals, self._fft_vals, axes=(2,3),
                threads=multiprocessing.cpu_count())
            self._ifft_plan = pyfftw.FFTW(
                self._fft_vals, self._vals, axes=(2,3),
                threads=multiprocessing.cpu_count(), direction="FFTW_BACKWARD")

            for i in range(0,len(self._wavelengths)):
                E_inputX = \
                    vn.vtk_to_numpy(point_data.GetArray(
                        "E_real_inputX_%sum" % self._wavelengths[i])) + \
                    1j*vn.vtk_to_numpy(point_data.GetArray(
                        "E_imag_inputX_%sum" % self._wavelengths[i]))
                E_inputY = \
                    vn.vtk_to_numpy(point_data.GetArray(
                        "E_real_inputY_%sum" % self._wavelengths[i])) + \
                    1j*vn.vtk_to_numpy(point_data.GetArray(
                        "E_imag_inputY_%sum" % self._wavelengths[i]))
                self._vals[i,[0,1],:,:] = \
                    E_inputX[:,[0,1]].transpose().reshape(2,dims[1],dims[0])
                self._vals[i,[2,3],:,:] = \
                    E_inputY[:,[0,1]].transpose().reshape(2,dims[1],dims[0])

            self._wavelengths = np.array(self._wavelengths).astype(float)
            (self._Nx, self._Ny) = (dims[0], dims[1])
            (self._dx, self._dy) = (spacings[0], spacings[1])
            (self._Lx, self._Ly) = (spacings[0]*(dims[0]-1), spacings[1]*(dims[1]-1))

        if len(kwargs)==3:
            if len(kwargs["mesh_dimensions"])!=2:
                raise Exception("mesh_dimensions should be an array-like object of size 2")
            if len(kwargs["mesh_lengths"])!=2:
                raise Exception("mesh_lengths should be an array-like object of size 2")

            self._wavelengths = np.array(kwargs["wavelengths"])

            dims = np.array(kwargs["mesh_dimensions"])
            lengths = np.array(kwargs["mesh_lengths"])

            (self._Nx, self._Ny) = tuple(int(dim) for dim in dims)
            (self._Lx, self._Ly) = tuple(lengths)
            (self._dx, self._dy) = tuple(lengths/np.maximum(np.ones(2),dims-1))

            self._vals = pyfftw.empty_aligned(
                (len(self._wavelengths),4,dims[1],dims[0]), dtype="complex128")
            self._fft_vals = pyfftw.empty_aligned(
                (len(self._wavelengths),4,dims[1],dims[0]), dtype="complex128")
            self._fft_plan = pyfftw.FFTW(
                self._vals, self._fft_vals, axes=(2,3),
                threads=multiprocessing.cpu_count())
            self._ifft_plan = pyfftw.FFTW(
                self._fft_vals, self._vals, axes=(2,3),
                threads=multiprocessing.cpu_count(), direction="FFTW_BACKWARD")
        
        IY, IX = np.meshgrid(range(0,self._Ny), range(0,self._Nx), indexing="ij")
        kx = -np.abs(2*np.pi/self._Lx*(IX-0.5*self._Nx)) + np.pi*self._Nx/self._Lx
        ky = -np.abs(2*np.pi/self._Ly*(IY-0.5*self._Ny)) + np.pi*self._Ny/self._Ly
        kr = np.tile(np.sqrt(kx**2+ky**2), (len(self._wavelengths),1,1)).flatten()
        k0 = np.tile(2*np.pi/self._wavelengths, (self._Nx,self._Ny,1)).transpose().flatten()

        filt = 1j*np.zeros(len(self._wavelengths)*self._Nx*self._Ny)
        filt[kr<k0] = np.exp(1j*np.sqrt(k0[kr<k0]**2-kr[kr<k0]**2))
        self._propag_filter = np.reshape(filt, (len(self._wavelengths),1,self._Ny,self._Nx))
        self._z = 0


    def copy(self):
        """Returns a hard copy of this OpticalFields object"""
        new_fields = OpticalFields(
            wavelengths = self._wavelengths,
            mesh_dimensions = (self._Nx, self._Ny),
            mesh_lengths = (self._Lx, self._Ly))
        new_fields.vals = self.vals # need to use the setter method for byte-aligned hard copy
        return new_fields


    @property
    def vals(self):
        """Numpy array for the optical fields values, of shape (N_wavelengths,4,Ny,Nx).

        If you want to initialize by hand the optical fields, the four components in the
        second dimension correspond to:
        
        * complex Ex field for an input polarisation//x
        * complex Ey field for an input polarisation//x
        * complex Ex field for an input polarisation//y
        * complex Ey field for an input polarisation//y

        """
        return self._vals


    @vals.setter
    def vals(self, fields_ndarray):
        if self._vals.shape==fields_ndarray.shape:
            self._vals[:] = fields_ndarray
        else:
            raise Exception("Wrong shape for the optical field ndarray")


    @property
    def z_focus(self):
        return self._z


    def propagate(self, new_z):
        """Propagate the optical fields to a new transverse plane at altitude ``new_z``
        by switching to Fourier space and using the exact propagator for the Helmholtz
        equation in free space."""
        filt = self._propag_filter**np.abs(new_z-self._z)
        filt = np.conj(filt) if new_z<self._z else filt

        self._fft_plan()
        self._fft_vals[:] *= filt
        self._ifft_plan()
        self._z = new_z


    def save_to_vti(self, filename):
        """Save the optical fields into a vti file.

        The ".vti" extension is automatically appended, no need to include it in the filename
        parameter (but in case you do only one extension will be added)
        """
        if filename[-4:]==".vti":
            path = filename
        else:
            path = filename+".vti"

        print("{ Saving optical fields to "+path+" }")
        vti_data = vtkImageData()
        vti_data.SetDimensions(self._Nx, self._Ny, 1)
        vti_data.SetOrigin(-self._Lx/2, -self._Ly/2, 0)
        vti_data.SetSpacing(self._dx, self._dy, 0)

        for i in range(0,len(self._wavelengths)):
            E_inputX = \
                self._vals[i,[0,1],:,:].reshape((2,self.get_n_vertices())).transpose()
            E_inputY = \
                self._vals[i,[2,3],:,:].reshape((2,self.get_n_vertices())).transpose()

            E_real_inputX = vn.numpy_to_vtk(np.real(E_inputX))
            E_real_inputX.SetName("E_real_inputX_%sum" % self._wavelengths[i])
            vti_data.GetPointData().AddArray(E_real_inputX)
            E_imag_inputX = vn.numpy_to_vtk(np.imag(E_inputX))
            E_imag_inputX.SetName("E_imag_inputX_%sum" % self._wavelengths[i])
            vti_data.GetPointData().AddArray(E_imag_inputX)

            E_real_inputY = vn.numpy_to_vtk(np.real(E_inputY))
            E_real_inputY.SetName("E_real_inputY_%sum" % self._wavelengths[i])
            vti_data.GetPointData().AddArray(E_real_inputY)
            E_imag_inputY = vn.numpy_to_vtk(np.imag(E_inputY))
            E_imag_inputY.SetName("E_imag_inputY_%sum" % self._wavelengths[i])
            vti_data.GetPointData().AddArray(E_imag_inputY)

        writer = vtkXMLImageDataWriter()
        writer.SetFileName(path)
        writer.SetInputData(vti_data)
        writer.Write()


    def get_pos(self, ix, iy):
        """Returns the position associated with the mesh indices (ix,iy)

        It is assumed that the mesh is centered on the origin (0,0).
        """
        return (ix*self._dx-self._Lx/2, iy*self._dy-self._Ly/2)


    def get_wavelengths(self):
        """Returns the wavelength array"""
        return self._wavelengths

    
    def get_mesh_dimensions(self):
        """Returns the dimensions (Nx,Ny) of the transverse mesh"""
        return (self._Nx, self._Ny)

    
    def get_mesh_lengths(self):
        """Returns the lengths (Lx,Ly) of the transverse mesh"""
        return (self._Lx, self._Ly)


    def get_mesh_spacings(self):
        """Returns the spacings (dx,dy,dz) of the transverse mesh"""
        return (self._dx, self._dy)


    def get_n_vertices(self):
        """Returns the number of vertices in the transverse mesh"""
        return self._Nx*self._Ny
