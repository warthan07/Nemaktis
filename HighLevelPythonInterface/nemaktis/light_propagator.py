import re
import os

import sys
sys.path.insert(0,"/home/gpoy/Nemaktis/BeamPropagationSolver/bin/")

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
    """The LightPropagator class allows to propagate optical fields through a LC sample as
    in a real microscope: a set of plane waves with different wavevectors and wavelengths
    are sent on the LC sample, and the associated transmitted optical fields (which can now
    longer be represented as plane waves due to diffraction) are calculated using one of the
    backend. 

    The actual set of wavelengths for the plane waves (choosen at construction) approximate
    the relevant part of the spectrum of the illumination light, whereas the set of
    wavevectors (also calculated at construction) are determined from the numerical aperture
    of the input condenser. The more open the condenser aperture is, the smoother the
    micrograph will look, since an open condenser aperture is associated with a wide range
    of angle for the wavectors of the incident plane waves. Conversely, an almost closed
    condenser aperture is associated with a single plane wave incident normally on the
    sample.

    Note that with the FieldViewer class, the transmitted optical fields calculated with
    this class can be projected on a visualisation screen through an objective of given
    numerical aperture. The numerical apertures of both the objective and condenser aperture
    can be set interactively in the FieldViewer class, whereas in this class we only
    specify the maximum value allowed for both quantities.

    The simulation and choice of backend is done by calling the method ``propagate_field``. 

    For each wavelength and wavevector of the incident plane wave, two simulations are done:
    one with a light source polarised along x, and one with a light source polarised along
    y. This allows us to fully caracterize the transmission of the LC sample and reconstruct
    any kind of optical micrograph.

    Parameters
    ----------
    material : :class:`~nemaktis.lc_material.LCMaterial` object
    wavelengths : array-like object
        An array containing all the wavelengths of the spectrum for the light source.
    max_NA_objective : float
        Sets the maximal numerical aperture for the microscope objective (you can
        dynamically adjust this quantity later on with a FieldViewer).
    max_NA_condenser : float
        Sets the maximal numerical aperture for the microscope condenser (you can
        dynamically adjust this quantity later on with a FieldViewer).
    N_radial_wavevectors : int
        Sets the number of wavevectors in the radial direction for the illumination plane
        waves. The total number of plane waves for each wavelength is 1+3*Nr*(Nr-1), where
        Nr correspond to the value of this parameter.
    """
    def __init__(self, *, material, wavelengths, max_NA_objective,
            max_NA_condenser = 0, N_radial_wavevectors = 1):
        if not isinstance(material, LCMaterial):
            raise TypeError("material should be a LCMaterial object")
        self._material = material
        self._wavelengths = list(wavelengths)
        self._max_NA_objective = max_NA_objective
        self._max_NA_condenser = max_NA_condenser
        self._N_radial_wavevectors = N_radial_wavevectors

        self._wavevectors = np.zeros(
            (1+3*N_radial_wavevectors*(N_radial_wavevectors-1),2))
        for ir in range(1,N_radial_wavevectors):
            beta = ir*self._max_NA_condenser/(N_radial_wavevectors-1)
            for iphi in range(0,6*ir):
                phi = iphi*np.pi/(3*ir)
                self._wavevectors[1+3*ir*(ir-1)+iphi,0] = beta*np.cos(phi)
                self._wavevectors[1+3*ir*(ir-1)+iphi,1] = beta*np.sin(phi)

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
        wavevectors = self._wavevectors.flatten().tolist()
        json_str = JSONEncoder().encode({
            "Algorithm settings": {
                "General": {
                    "LC field type":               "Director",
                    "Results folder name":         "" },
                "Beam propagation": {
                    "N Woodbury steps":            2,
                    "Number of substeps per slab": 1 }},
            "Physics settings": {
                "Initial conditions": {
                    "Beam profile":     "UniformBeam",
                    "LC field file":    "",
                    "Mesh dimensions":  dims,
                    "Mesh spacings":    spacings,
                    "Basis convention": "XYZ" },
                "Coefficients": {
                    "no":               str(self._material.no),
                    "ne":               str(self._material.ne),
                    "nhost":            str(self._material.nhost),
                    "nin":              str(self._material.nin),
                    "Wavelengths":      self._wavelengths,
                    "Wavevectors":      wavevectors}},
            "Postprocessor settings": {
                "Bulk output": {
                    "Activate":	        bulk_filename is not None,
                    "Base name":        bulk_filename if bulk_filename is not None else ""},
                "Screen output": {
                    "Activate":	                          True,
                    "Base name":                          "",
                    "Isotropic layer thicknesses":        self._material.iso_layer_thicknesses,
                    "Isotropic layer refractive indices": self._material.iso_layer_indices,
                    "Focalisation z-shift":               0,
                    "Numerical aperture":                 self._max_NA_objective }}})

        director_vals = director_field.vals.ravel()
        if director_field.mask_type is not None:
            mask_vals = director_field.mask_vals.ravel()

        N_fields_vals = \
            self._wavevectors.shape[0]*len(self._wavelengths)*4*dims[0]*dims[1]
        if director_field.mask_type is None:
            data_out = bpm.run_backend_without_mask(
                json_str, director_vals, N_fields_vals)
        else:
            data_out = bpm.run_backend_with_mask(
                json_str, director_vals, mask_vals, N_fields_vals)
        print("")

        output_fields = OpticalFields(
            wavelengths = self._wavelengths,
            max_NA_objective = self._max_NA_objective,
            max_NA_condenser = self._max_NA_condenser,
            N_radial_wavevectors = self._N_radial_wavevectors,
            mesh_lengths = (spacings[0]*(dims[0]-1), spacings[1]*(dims[1]-1)),
            mesh_dimensions = (dims[0], dims[1]))

        Nl = len(self._wavelengths)
        Nq = len(self._wavevectors)
        output_fields.vals = data_out.reshape((Nl, Nq, 4, dims[1], dims[0]))/np.sqrt(2)
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
                print("Warning: dtmm does not support dispersive index; " +
                      "Using ne(0.6µm) instead")
            l = 0.6
            ne = eval(self._material.ne.replace("lambda","l").replace("^","**"))
        else:
            ne = self._material.ne
        if isinstance(self._material.no, str):
            if "lambda" in self._material.no:
                print("Warning: dtmm does not support dispersive index; " +
                      "Using no(0.6µm) instead")
            l = 0.6
            no = eval(self._material.no.replace("lambda","l").replace("^","**"))
        else:
            no = self._material.no
        if isinstance(self._material.nhost, str):
            if "lambda" in self._material.nhost:
                print("Warning: dtmm does not support dispersive index; " +
                      "Using nhost(0.6µm) instead")
            l = 0.6
            nhost = eval(self._material.nhost.replace("lambda","l").replace("^","**"))
        else:
            nhost = self._material.nhost
        if isinstance(self._material.nhost, str):
            if "lambda" in self._material.nin:
                print("Warning: dtmm does not support dispersive index; " +
                      "Using nin(0.6µm) instead")
            l = 0.6
            nin = eval(self._material.nin.replace("lambda","l").replace("^","**"))
        else:
            nin = self._material.nin

        if len(self._material.iso_layer_indices)!=0:
            print("Warning: specified isotropic layers will be ignored since this feature is "+
                  "not yet supported in dtmm.")
        print("")
        
        if director_field.mask_vals is not None:
            mask_vals = director_field.mask_vals>=0
        else:
            mask_vals = None
        optical_data = dtmm.director2data(
           director_field.vals, mask = mask_vals,
           no = no, ne = ne, nhost = nhost,
           thickness = spacings[2]/spacings[0]*np.ones(dims[2]))

        wavelengths = 1000*np.array(self._wavelengths)
        beta = np.zeros((len(self._wavevectors),))
        phi = np.zeros((len(self._wavevectors),))
        intensity = np.ones((len(self._wavevectors),))
        for ir in range(1,self._N_radial_wavevectors):
            for iphi in range(0,6*ir):
                beta[1+3*ir*(ir-1)+iphi] = ir*self._max_NA_condenser/(self._N_radial_wavevectors-1)
                phi[1+3*ir*(ir-1)+iphi] = iphi*np.pi/3

        field_data_in = dtmm.illumination_data(
            (dims[1],dims[0]), wavelengths, pixelsize=1000*spacings[0], n=nin,
            beta=beta, phi=phi, intensity=intensity)
        field_data_out = dtmm.transfer_field(
            field_data_in, optical_data, nin=nin,
            betamax=self._max_NA_objective, diffraction=diffraction,
            ret_bulk=bulk_filename is not None)[0]
        print("")

        if bulk_filename is not None:
            print("{ Saving optical fields to "+bulk_filename+".vti }")
            lengths = director_field.get_mesh_lengths()

            vti_data = vtkImageData()
            vti_data.SetDimensions(dims[0], dims[1], dims[2]+1)
            vti_data.SetOrigin(-lengths[0]/2, -lengths[1]/2, -lengths[2]/2)
            vti_data.SetSpacing(spacings[0], spacings[1], spacings[2]*(dims[2]-1)/dims[2])

            wavelengths_data = vn.numpy_to_vtk(self._wavelengths)
            wavelengths_data.SetName("lambda")
            vti_data.GetFieldData().AddArray(wavelengths_data)

            qx_data = vn.numpy_to_vtk(self._wavevectors[:,0])
            qx_data.SetName("qx")
            vti_data.GetFieldData().AddArray(qx_data)

            qy_data = vn.numpy_to_vtk(self._wavevectors[:,1])
            qy_data.SetName("qy")
            vti_data.GetFieldData().AddArray(qy_data)

            Np = dims[0]*dims[1]*(dims[2]+1)
            for wave_idx in range(0,len(self._wavelengths)):
                for q_idx in range(0,len(self._wavevectors)):
                    E_inputX = field_data_out[:-1,q_idx,0,wave_idx,[0,2],:,:].transpose(
                        (1,0,2,3)).reshape((2,Np)).transpose()
                    E_inputY = field_data_out[:-1,q_idx,1,wave_idx,[0,2],:,:].transpose(
                        (1,0,2,3)).reshape((2,Np)).transpose()

                    E_real_inputX = vn.numpy_to_vtk(np.real(E_inputX))
                    E_real_inputX.SetName("E_real_inputX_%d_%d" % (wave_idx,q_idx))
                    vti_data.GetPointData().AddArray(E_real_inputX)
                    E_imag_inputX = vn.numpy_to_vtk(np.imag(E_inputX))
                    E_imag_inputX.SetName("E_imag_inputX_%d_%d" % (wave_idx,q_idx))
                    vti_data.GetPointData().AddArray(E_imag_inputX)

                    E_real_inputY = vn.numpy_to_vtk(np.real(E_inputY))
                    E_real_inputY.SetName("E_real_inputY_%d_%d" % (wave_idx,q_idx))
                    vti_data.GetPointData().AddArray(E_real_inputY)
                    E_imag_inputY = vn.numpy_to_vtk(np.imag(E_inputY))
                    E_imag_inputY.SetName("E_imag_inputY_%d_%d" % (wave_idx,q_idx))
                    vti_data.GetPointData().AddArray(E_imag_inputY)

            writer = vtkXMLImageDataWriter()
            writer.SetFileName(bulk_filename+".vti")
            writer.SetInputData(vti_data)
            writer.Write()

            # We only keep the last slice to compute micrographs
            field_data_out = field_data_out[-1,:,:,:,:,:,:]

        output_fields = OpticalFields(
            wavelengths = self._wavelengths,
            max_NA_objective = self._max_NA_objective,
            max_NA_condenser = self._max_NA_condenser,
            N_radial_wavevectors = self._N_radial_wavevectors,
            mesh_lengths = (spacings[0]*(dims[0]-1), spacings[1]*(dims[1]-1)),
            mesh_dimensions = (dims[0], dims[1]))

        Nl = len(self._wavelengths)
        Nq = len(self._wavevectors)
        output_fields.vals = field_data_out[:,:,:,[0,2],:,:].transpose(
            (2,0,1,3,4,5)).reshape((Nl,Nq,4,dims[1],dims[0]))
        return output_fields


class OpticalFields:
    """The OpticalFields object stores the mesh information of the transverse mesh (plane mesh
    orthogonal to the z-direction, default altitude of 0) and the optical fields values on
    this mesh.  Since this python package is mainly used to reconstruct micrographs, we only
    store internally the complex horizontal electric field for two simulation: one with a
    light source polarised along ``x``, and the other with a light source polarised along
    ``y``.  In case multiple wavelengths/wavectors were used in the simulation, we store these
    quantities separately for each wavelength/wavevector.

    This class is initialised either manually or with a path to a vti file containing
    previously calculated optical fields and mesh details.
    
    In the first version of this constructor:

    .. code-block:: python
    
        optical_fields = OpticalFields(
            wavelengths=[l0,l1,...,lN], max_NA_objective=NA_o,
            max_NA_condenser=NA_c, N_radial_wavevectors=Nr, 
            mesh_lengths=(Lx,Ly), mesh_dimensions=(Nx,Ny))

    the actual values of the transverse fields needs to be provided later using the raw
    setter method fields_vals (shape (N_wavelengths,N_wavevectors,4,Ny,Nx), with
    N_wavevectors=3*Nr*(Nr-1)+1).

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
            field_data = reader.GetOutput().GetFieldData()
            dims = np.array(reader.GetOutput().GetDimensions())
            spacings = np.array(reader.GetOutput().GetSpacing())

            if dims[2]!=1:
                raise Exception("The specified vti file should include 2D data")

            if not field_data.HasArray("lambda"):
                raise Exception(
                    "VTI file is missing the field array \"lambda\" for the wavelengths")
            self._wavelengths = vn.vtk_to_numpy(field_data.GetArray("lambda"))

            if not field_data.HasArray("qx"):
                raise Exception(
                    "VTI file is missing the field array \"qx\" for the wavevectors")
            if not field_data.HasArray("qy"):
                raise Exception(
                    "VTI file is missing the field array \"qy\" for the wavevectors")
            self._wavevectors = np.stack(
                (vn.vtk_to_numpy(field_data.GetArray("qx")),
                vn.vtk_to_numpy(field_data.GetArray("qy"))), axis=-1)
            Nq = len(self._wavevectors)
            Nr = int(np.round((1+np.sqrt((4*Nq-1)/3))/2))
            if Nq!=1+3*Nr*(Nr-1):
                raise Exception(
                    "VTI file contain the wrong number of wavevectors")
            else:
                self._N_radial_wavevectors = Nr
            self._max_NA_condenser = np.sqrt(np.sum(self._wavevectors**2,axis=1))[-1]
            for qr_idx in range(0,Nr):
                q_idx_start = 1+3*qr_idx*(qr_idx-1) if qr_idx>0 else 0
                q_idx_end = 1+3*qr_idx*(qr_idx+1)
                q_norms = np.sqrt(np.sum(self._wavevectors[q_idx_start:q_idx_end,:]**2,axis=1))
                if np.max(np.abs(q_norms-qr_idx*self._max_NA_condenser/(Nr-1)))>1e-8:
                    raise Exception(
                        "Incompatible wavevector mesh inside the VTI file")

            if not field_data.HasArray("max_NA_objective"):
                raise Exception(
                    "VTI file is missing the field scalar \"max_NA_objective\"")
            self._max_NA_objective = \
                vn.vtk_to_numpy(field_data.GetArray("max_NA_objective"))[0]

            Nl = len(self._wavelengths)
            self._vals = pyfftw.empty_aligned(
                (Nl,Nq,4,dims[1],dims[0]), dtype="complex128")
            self._fft_vals = pyfftw.empty_aligned(
                (Nl,Nq,4,dims[1],dims[0]), dtype="complex128")
            self._focused_vals = pyfftw.empty_aligned(
                (Nl,Nq,4,dims[1],dims[0]), dtype="complex128")
            self._fft_plan = pyfftw.FFTW(
                self._vals, self._fft_vals, axes=(3,4),
                threads=multiprocessing.cpu_count())
            self._ifft_plan = pyfftw.FFTW(
                self._fft_vals, self._focused_vals, axes=(3,4),
                threads=multiprocessing.cpu_count(), direction="FFTW_BACKWARD")

            for wave_idx in range(0,Nl):
                for q_idx in range(0,Nq):
                    suffixes = ["real_inputX", "imag_inputX", "real_inputY", "imag_inputY"]
                    for suffix in suffixes:
                        array_name = "E_"+suffix+"_%d_%d" % (wave_idx,q_idx)
                        if not point_data.HasArray(array_name):
                            raise Exception(
                                "Missing array \""+array_name+"\" in VTI file")
                    E_inputX = \
                        vn.vtk_to_numpy(point_data.GetArray(
                            "E_real_inputX_%d_%d" % (wave_idx,q_idx))) + \
                        1j*vn.vtk_to_numpy(point_data.GetArray(
                            "E_imag_inputX_%d_%d" % (wave_idx,q_idx)))
                    E_inputY = \
                        vn.vtk_to_numpy(point_data.GetArray(
                            "E_real_inputY_%d_%d" % (wave_idx,q_idx))) + \
                        1j*vn.vtk_to_numpy(point_data.GetArray(
                            "E_imag_inputY_%d_%d" % (wave_idx,q_idx)))
                    self._vals[wave_idx,q_idx,[0,1],:,:] = \
                        E_inputX[:,[0,1]].transpose().reshape(2,dims[1],dims[0])
                    self._vals[wave_idx,q_idx,[2,3],:,:] = \
                        E_inputY[:,[0,1]].transpose().reshape(2,dims[1],dims[0])
    
            (self._Nx, self._Ny) = (dims[0], dims[1])
            (self._dx, self._dy) = (spacings[0], spacings[1])
            (self._Lx, self._Ly) = (spacings[0]*(dims[0]-1), spacings[1]*(dims[1]-1))

        else:
            if len(kwargs["mesh_dimensions"])!=2:
                raise Exception("mesh_dimensions should be an array-like object of size 2")
            if len(kwargs["mesh_lengths"])!=2:
                raise Exception("mesh_lengths should be an array-like object of size 2")

            self._wavelengths = np.array(kwargs["wavelengths"])
            self._max_NA_objective = kwargs["max_NA_objective"]
            self._max_NA_condenser = kwargs["max_NA_condenser"]
            self._N_radial_wavevectors = kwargs["N_radial_wavevectors"]

            Nr = self._N_radial_wavevectors
            self._wavevectors = np.zeros((1+3*Nr*(Nr-1),2))
            for ir in range(1,Nr):
                beta = ir*self._max_NA_condenser/(Nr-1)
                for iphi in range(0,6*ir):
                    phi = iphi*np.pi/(3*ir)
                    self._wavevectors[1+3*ir*(ir-1)+iphi,0] = beta*np.cos(phi)
                    self._wavevectors[1+3*ir*(ir-1)+iphi,1] = beta*np.sin(phi)

            dims = np.array(kwargs["mesh_dimensions"])
            lengths = np.array(kwargs["mesh_lengths"])

            (self._Nx, self._Ny) = tuple(int(dim) for dim in dims)
            (self._Lx, self._Ly) = tuple(lengths)
            (self._dx, self._dy) = tuple(lengths/np.maximum(np.ones(2),dims-1))

            Nl = len(self._wavelengths)
            Nq = len(self._wavevectors)
            self._vals = pyfftw.empty_aligned(
                (Nl,Nq,4,dims[1],dims[0]), dtype="complex128")
            self._fft_vals = pyfftw.empty_aligned(
                (Nl,Nq,4,dims[1],dims[0]), dtype="complex128")
            self._focused_vals = pyfftw.empty_aligned(
                (Nl,Nq,4,dims[1],dims[0]), dtype="complex128")
            self._fft_plan = pyfftw.FFTW(
                self._vals, self._fft_vals, axes=(3,4),
                threads=multiprocessing.cpu_count())
            self._ifft_plan = pyfftw.FFTW(
                self._fft_vals, self._focused_vals, axes=(3,4),
                threads=multiprocessing.cpu_count(), direction="FFTW_BACKWARD")

        if self._N_radial_wavevectors>1 and self._max_NA_condenser>0:
            self._delta_qr = self._max_NA_condenser/(self._N_radial_wavevectors-1)
        else:
            self._delta_qr = 1


        IY, IX = np.meshgrid(range(0,self._Ny), range(0,self._Nx), indexing="ij")
        kx = -np.abs(2*np.pi/self._Lx*(IX-0.5*self._Nx)) + np.pi*self._Nx/self._Lx
        ky = -np.abs(2*np.pi/self._Ly*(IY-0.5*self._Ny)) + np.pi*self._Ny/self._Ly

        k0 = np.tile(2*np.pi/self._wavelengths, (self._Nx,self._Ny,1,Nq,1)).transpose()
        px = k0*self._wavevectors[:,0][np.newaxis,:,np.newaxis,np.newaxis,np.newaxis]
        py = k0*self._wavevectors[:,1][np.newaxis,:,np.newaxis,np.newaxis,np.newaxis]

        kSqr = (kx[np.newaxis,np.newaxis,:,:]+px)**2+(ky[np.newaxis,np.newaxis,:,:]+py)**2
        mask = kSqr.flatten()<k0.flatten()**2

        filt = 1j*np.zeros(Nl*Nq*self._Nx*self._Ny)
        filt[mask] = np.exp(1j*np.sqrt(k0.flatten()[mask]**2-kSqr.flatten()[mask]))
        self._objective_filter = np.reshape(filt, (Nl,Nq,1,self._Ny,self._Nx))
        self._objective_mask = (kSqr<(k0*self._max_NA_objective)**2).astype(float)
        self._kSqr = kSqr
        self._k0 = k0
        self._z = 0


    def copy(self):
        """Returns a hard copy of this OpticalFields object"""
        new_fields = OpticalFields(
            wavelengths = self._wavelengths,
            max_NA_condenser = self._max_NA_condenser,
            N_radial_wavevectors = self._N_radial_wavevectors,
            mesh_dimensions = (self._Nx, self._Ny),
            mesh_lengths = (self._Lx, self._Ly))
        new_fields.vals = self.vals # need to use the setter method for byte-aligned hard copy
        return new_fields


    @property
    def focused_vals(self):
        """Numpy array for the optical fields values after focalisation by the microscope
        objective, of shape (N_wavelengths,N_wavevectors,4,Ny,Nx).
        """
        return self._focused_vals


    @property
    def vals(self):
        """Numpy array for the optical fields values, of shape 
        (N_wavelengths,N_wavevectors,4,Ny,Nx).

        If you want to initialize by hand the optical fields, the four components in the
        third dimension correspond to:
        
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


    def update_NA_objective(self, new_NA):
        NA = max(0,min(self._max_NA_objective,new_NA))
        self._objective_mask = (self._kSqr<(self._k0*NA)**2).astype(float)


    def focus_fields(self, z_focus=None):
        """Propagate the optical fields through the objective lens to the screen conjugate
        to the focusing plane (whose altitude inside the sample is set with the parameter
        z_focus)."""
        if z_focus is not None:
            filt = self._objective_mask*self._objective_filter**np.abs(z_focus)
            self._z = z_focus
        else:
            z_focus = self._z
            filt = self._objective_mask*self._objective_filter**np.abs(z_focus)

        self._fft_plan()
        self._fft_vals *= np.conj(filt) if z_focus<0 else filt
        self._ifft_plan()


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

        Np = self.get_n_vertices()
        for wave_idx in range(0,len(self._wavelengths)):
            for q_idx in range(0,len(self._wavevectors)):
                E_inputX = self._vals[wave_idx,q_idx,[0,1],:,:].reshape((2,Np)).transpose()
                E_inputY = self._vals[wave_idx,q_idx,[2,3],:,:].reshape((2,Np)).transpose()

                E_real_inputX = vn.numpy_to_vtk(np.real(E_inputX))
                E_real_inputX.SetName("E_real_inputX_%d_%d" % (wave_idx,q_idx))
                vti_data.GetPointData().AddArray(E_real_inputX)
                E_imag_inputX = vn.numpy_to_vtk(np.imag(E_inputX))
                E_imag_inputX.SetName("E_imag_inputX_%d_%d" % (wave_idx,q_idx))
                vti_data.GetPointData().AddArray(E_imag_inputX)

                E_real_inputY = vn.numpy_to_vtk(np.real(E_inputY))
                E_real_inputY.SetName("E_real_inputY_%d_%d" % (wave_idx,q_idx))
                vti_data.GetPointData().AddArray(E_real_inputY)
                E_imag_inputY = vn.numpy_to_vtk(np.imag(E_inputY))
                E_imag_inputY.SetName("E_imag_inputY_%d_%d" % (wave_idx,q_idx))
                vti_data.GetPointData().AddArray(E_imag_inputY)

        wavelengths_data = vn.numpy_to_vtk(self._wavelengths)
        wavelengths_data.SetName("lambda")
        vti_data.GetFieldData().AddArray(wavelengths_data)

        qx_data = vn.numpy_to_vtk(self._wavevectors[:,0])
        qx_data.SetName("qx")
        vti_data.GetFieldData().AddArray(qx_data)

        qy_data = vn.numpy_to_vtk(self._wavevectors[:,1])
        qy_data.SetName("qy")
        vti_data.GetFieldData().AddArray(qy_data)

        NA_data = vn.numpy_to_vtk(np.array([self._max_NA_objective]))
        NA_data.SetName("max_NA_objective")
        vti_data.GetFieldData().AddArray(NA_data)

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


    def get_wavevectors(self):
        """Returns the wavevectors array"""
        return self._wavevectors


    def get_qr_index(self, NA_condenser):
        """For internal use.

        Allows to build sub-range of wavevector index for a given numerical aperture of the
        condenser, which must be smaller than the internal maximal numerical aperture set at
        construction.
        """
        return int(np.ceil(NA_condenser/self._delta_qr))

    def get_delta_qr(self):
        """For internal use.
        
        Allows to build integration rule with respect to the wavectors.
        """
        return self._delta_qr


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
