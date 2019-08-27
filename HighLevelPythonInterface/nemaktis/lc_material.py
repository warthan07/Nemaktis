import numpy as np
import os

from vtk import vtkImageData, vtkXMLImageDataReader, vtkXMLImageDataWriter
from vtk.util import numpy_support as vn

from warnings import simplefilter
simplefilter(action='ignore', category=FutureWarning)


class LCMaterial(object):
    """A class containing the director field data, simulation mesh, and
    physics constants.

    Parameters
    ----------
    director_field : :class:`~nemaktis.lc_material.DirectorField` object
    ne : float or math string depending on the wavelength variable "lambda" (µm)
        The extraordinary refractive index associated with the LC
        material.
    no : float or math string depending on the wavelength variable "lambda" (µm)
        The ordinary refractive index associated with the LC material.
    nhost : optional, float or math string depending on the wavelength variable "lambda" (µm)
        The refractive index associated with an eventual host fluid in which the LC domain
        is embedded (see DirectorField.set_mask).
    nin : optional, float or math string depending on the wavelength variable "lambda" (µm)
        The refractive index associated with the input medium below the LC layer. 
        A default value of 1 is assumed.
    """
    iso_layer_indices = [] 
    iso_layer_thicknesses = [] 

    def __init__(self, *, director_field, ne, no, nhost = 1, nin = 1):
        self.director_field = director_field
        self.ne = ne
        self.no = no
        self.nhost = nhost
        self.nin = nin

    def add_isotropic_layer(self, *, nlayer, thickness):
        """Add an isotropic layer above the sample.  Light is assumed to propagate in the
        z-direction, and will cross first the LC material, and then the isotropic layers
        specified with this function.

        Parameters
        ----------
        nlayer : float
            Refractive index of the new isotropic layer
        thickness : float
            Thickness (µm) of the new isotropic layer
        """
        self.iso_layer_indices.append(nlayer)
        self.iso_layer_thicknesses.append(thickness)



class DirectorField:
    """The DirectorField class stores the director field and mesh
    informations. It is initialized given either the lengths and
    dimensions of the associated 3D mesh or a path to a vti file
    containing the director field and mesh details.

    In the first version of this constructor:

    .. code-block:: python
    
        nfield = DirectorField(
            mesh_lengths=(Lx,Ly,Lz), mesh_dimensions=(Nx,Ny,Nz))

    the actual values of the director field needs to be provided later
    using the "init_from_funcs" method or via the setter method "vals"
    (numpy array of shape (Nz,Ny,Nx,3)). The mesh lengths needs to be
    specified in micrometer.

    In the second version of this constructor:

    .. code-block:: python
    
        nfield = DirectorField(vti_file="path to vti file")

    the values of the director field and the details of the mesh are
    automatically assigned from the vti file.
    """

    _mask_type = None
    _mask_formula = None
    _mask_vals = None

    def __init__(self, **kwargs):
        if len(kwargs)==1:
            if not os.path.isfile(kwargs["vti_file"]):
                raise Exception("VTI file does not exists")

            print("{ Initializing director field from "+kwargs["vti_file"]+" }")

            reader = vtkXMLImageDataReader()
            reader.SetFileName(kwargs["vti_file"])
            reader.Update()

            vti_data = reader.GetOutput()
            dims = np.array(vti_data.GetDimensions())
            spacings = np.array(vti_data.GetSpacing())

            if not vti_data.GetPointData().HasArray("n"):
                raise Exception("Missing vtk array 'n' in the given vti file")
            nfield = vn.vtk_to_numpy(vti_data.GetPointData().GetArray("n"))

            self._vals =  nfield.reshape((dims[2], dims[1], dims[0], 3))
            (self._Nx, self._Ny, self._Nz) = tuple(int(dim) for dim in dims)
            (self._dx, self._dy, self._dz) = tuple(spacings)
            (self._Lx, self._Ly, self._Lz) = tuple(spacings*(dims-1))

            if vti_data.GetPointData().HasArray("lc_mask"):
                mask_vals = vn.vtk_to_numpy(vti_data.GetPointData().GetArray("lc_mask"))
                self._mask_type = "raw"
                self._mask_vals = mask_vals.reshape((dims[2],dims[1],dims[0]))


        if len(kwargs)==2:
            if len(kwargs["mesh_lengths"])!=3:
                raise Exception("mesh_lengths should be an array-like object of size 3")
            if len(kwargs["mesh_dimensions"])!=3:
                raise Exception("mesh_dimensions should be an array-like object of size 3")

            print("{ Initializing empty director field }")

            dims = np.array(kwargs["mesh_dimensions"])
            lengths = np.array(kwargs["mesh_lengths"])

            (self._Nx, self._Ny, self._Nz) = tuple(int(dim) for dim in dims)
            (self._Lx, self._Ly, self._Lz) = tuple(lengths)
            (self._dx, self._dy, self._dz) = tuple(lengths/np.maximum(np.ones(3),dims-1))
            self._vals = np.zeros((self._Nz,self._Ny,self._Nx,3))


    def set_mask(self, *,  mask_type, mask_formula = None, mask_ndarray = None):
        """Set a boolean mask for the LC domain. This method allows to specifify complex shape
        for the LC domain inside the regular cartesian mesh specified at construction. 

        Three possible ways of initializing the mask are possible. If you simply want to
        specify a spherical domain for a droplet centered on the mesh and of diameter equal to
        the mesh length along z, call:

        .. code-block:: python
            
            nfield.set_mask(mask_type="droplet")

        You can also use a string formula depending on the space variables ``x``, ``y`` and
        ``z`` and which must evaluates to True if the associated point is inside the LC
        domain, else False:

        .. code-block:: python
            
            nfield.set_mask(mask_type="formula", mask_formula="your formula")

        Finally, you can directly gives a boolean numpy array of shape (Nz,Ny,Nx), where each
        value in this array must be True if the associated mesh point is inside the LC domain,
        else False:
        
        .. code-block:: python
            
            nfield.set_mask(mask_type="raw", mask_ndarray=your_mask_array)

        """
        if mask_type=="droplet":
            self._mask_type = mask_type
            self._mask_formula = str((self._Lz/2)**2)+">x**2+y**2+z**2"
        elif mask_type=="formula":
            self._mask_type = mask_type
            self._mask_formula = mask_formula.replace("^","**")
        elif mask_type=="raw":
            if(mask_ndarray.shape!=(self._Nz,self._Ny,self._Nx)):
                raise Exception("Wrong shape for the mask ndarray")
            self._mask_type = mask_type
            self._mask_vals = mask_ndarray
        else:
            raise Exception("Unrecognised mask type: should be 'droplet', 'formula' or 'raw'")


    def delete_mask(self):
        """Delete the current LC mask."""
        self._mask_type = None
        self._mask_formula = None
        self._mask_vals = None


    @property
    def mask_type(self):
        """Returns the mask type: "droplet", "formula" or "raw"."""
        return self._mask_type


    @property
    def mask_formula(self):
        """Returns the LC mask formula if it was set, else returns None."""
        return self._mask_formula


    @property
    def mask_vals(self):
        """Returns the LC mask boolean array."""
        if self._mask_type=="raw":
            return self._mask_vals
        elif self._mask_type=="formula" or self._mask_type=="droplet":
            mask_vals = np.ndarray((self._Nz,self._Ny,self._Nx))
            for iz in range(0,self._Nz):
                z = self._dz*(iz-0.5*(self._Nz-1))
                for iy in range(0,self._Ny):
                    y = self._dy*(iy-0.5*(self._Ny-1))
                    for ix in range(0,self._Nx):
                        x = self._dx*(ix-0.5*(self._Nx-1))
                        mask_vals[iz,iy,ix] = eval(self._mask_formula)
            return mask_vals
        else:
            return None

    
    def extend(self, scale_x, scale_y):
        """
        Extend the computational mesh in the ``xy`` plane by padding new points near the ``x``
        and ``y`` boundaries. These new points are initialized with the edge value of the
        director field on the ``x`` and ``y`` boundaries.

        Parameters
        ----------
        scale_x : float
            The mesh length in the x-direction will be scaled by this factor.
        scale_y : float
            The mesh length in the y-direction will be scaled by this factor.

        """
        print("{ Extending computational mesh }")

        Nx_pad = int(self._Nx*(scale_x-1)/2)
        Ny_pad = int(self._Ny*(scale_y-1)/2)
        self._vals = np.pad(
            self._vals, ((0,0),(Ny_pad,Ny_pad),(Nx_pad,Nx_pad),(0,0)), "edge")
        if self._mask_vals is not None:
            raise Exception("You should always call set_mask after extend")

        self._Nx = self._Nx+2*Nx_pad
        self._Ny = self._Ny+2*Ny_pad

        self._Lx = self._dx*(self._Nx-1)
        self._Ly = self._dy*(self._Ny-1)


    def rotate_90deg(self, axis):
        """
        Rotate the director field by 90 degrees around the specified axis.
        
        Parameters
        ----------
        axis : str
            Axis around which to perform the rotation. Need to be under the form '[s]A' where
            the optional parameter 's'='+' or '-' decribes the sign of rotation and 'A'='x',
            'y' or 'z' defines the rotation axis.
        """
        if axis[-1]=="x":
            ax = 0
        elif axis[-1]=="y":
            ax = 1
        elif axis[-1]=="z":
            ax = 2
        else:
            raise Exception("Could not parse axis.")

        if len(axis)==1 or (len(axis)==2 and axis[0]=="+"):
            flip_ax = (ax+1)%3
        elif len(axis)==2 and axis[0]=="-":
            flip_ax = (ax+2)%3
        else:
            raise Exception("Could not parse axis.")

        print("{ 90 deg rotation around axis "+axis[-1]+" }")

        # mapping for the axis order in the raw data array (axis 3 is component axis)
        ax_map = [2, 1, 0, 3]

        # permutation of axes and components equivalent to the rotation operation
        permut_ax = np.arange(0,4)
        permut_ax[ax_map[(ax+1)%3]] = ax_map[(ax+2)%3]
        permut_ax[ax_map[(ax+2)%3]] = ax_map[(ax+1)%3]

        permut_comp = np.arange(0,3)
        permut_comp[(ax+1)%3] = (ax+2)%3
        permut_comp[(ax+2)%3] = (ax+1)%3

        # We apply the permutations
        self._vals = np.flip(
            np.transpose(self.vals, tuple(permut_ax)),
            axis = ax_map[flip_ax])[:,:,:,permut_comp]
        self._vals[:,:,:,flip_ax] = -self._vals[:,:,:,flip_ax]

        (self._Nx, self._Ny, self._Nz) = tuple(
            self.get_mesh_dimensions()[i] for i in permut_comp[0:3])
        (self._Lx, self._Ly, self._Lz) = tuple(
            self.get_mesh_lengths()[i] for i in permut_comp[0:3])
        (self._dx, self._dy, self._dz) = tuple(
            self.get_mesh_spacings()[i] for i in permut_comp[0:3])


    def rotate_180deg(self, axis):
        """
        Rotate the director field by 180 degrees around the specified axis.
        
        Parameters
        ----------
        axis : str
            Axis around which to perform the rotation. Need to be under the form 'A' where
            'A'='x', 'y' or 'z' defines the rotation axis.
        """
        if axis[0]=="x":
            ax = 0
        elif axis[0]=="y":
            ax = 1
        elif axis[0]=="z":
            ax = 2
        else:
            raise Exception("Could not parse axis.")

        print("{ 180 deg rotation around axis "+axis[-1]+" }")

        # mapping for the axis order in the raw data array (axis 3 is component axis)
        ax_map = np.array([2, 1, 0, 3])

        # Axes to be flipped
        flip_axes = [(ax+1)%3, (ax+2)%3]

        # We apply the rotation
        self._vals = np.flip(self.vals, axis = ax_map[flip_axes])
        self._vals[:,:,:,flip_axes] = -self._vals[:,:,:,flip_axes]


    @property
    def vals(self):
        """Numpy array for the director values, of shape (Nz,Ny,Nx,3)."""
        return self._vals


    @vals.setter
    def vals(self, director_ndarray):
        if self._vals.shape==director_ndarray.shape:
            self._vals = director_ndarray
        else:
            raise Exception("Wrong shape for the director field ndarray")


    def init_from_funcs(self, nx_func, ny_func, nz_func):
        """Initialize the director field from three functions for each of its component. The
        functions must depends on the space variables ``x``, ``y`` and ``z``. We recall that
        the mesh is centered on the origin.

        If the given function are numpy-vectorizable, this function should be pretty fast. If
        not, a warning will be printed and the given function will be vectorized with the
        numpy method ``vectorize`` (in which case you should expect a much slower execution
        time).

        """

        print("{ Calculating director values from user functions }")
        zz, yy, xx = np.meshgrid(np.linspace(-self._Lz/2, self._Lz/2, self._Nz),
                                 np.linspace(-self._Ly/2, self._Ly/2, self._Ny),
                                 np.linspace(-self._Lx/2, self._Lx/2, self._Nx),
                                 indexing="ij")

        # We verify if the user functions are vectorizable
        dummy_arr = np.ones((2,2,2))
        try:
            nx_func(dummy_arr,dummy_arr,dummy_arr)
        except:
            print("\tnx_func is not vectorized, using a non-optimized version instead.")
            nx_func = np.vectorize(nx_func)
        try:
            ny_func(dummy_arr,dummy_arr,dummy_arr)
        except:
            print("\tny_func is not vectorized, using a non-optimized version instead.")
            ny_func = np.vectorize(ny_func)
        try:
            nz_func(dummy_arr,dummy_arr,dummy_arr)
        except:
            print("\tnz_func is not vectorized, using a non-optimized version instead.")
            nz_func = np.vectorize(nz_func)

        self._vals = np.concatenate((np.expand_dims(nx_func(xx, yy, zz), axis=3),
                                       np.expand_dims(ny_func(xx, yy, zz), axis=3),
                                       np.expand_dims(nz_func(xx, yy, zz), axis=3)), 3)

    def normalize(self):
        """Normalize the director field values to 1."""

        print("{ Normalizing director values to 1 }")
        norms = np.sqrt(np.sum(self._vals**2, axis=3, keepdims=True))
        norms[norms==0] = 1
        self._vals = self._vals / np.tile(norms, (1,1,1,3))
    

    def save_to_vti(self, filename):
        """Save the director field into a vti file.

        The ".vti" extension is automatically appended, no need to include it in the filename
        parameter (but in case you do only one extension will be added)
        """ 
        if filename[-4:]==".vti":
            path = filename
        else:
            path = filename+".vti"

        print("{ Saving director field to "+path+" }")
        vti_data = vtkImageData()
        vti_data.SetDimensions(self._Nx, self._Ny, self._Nz)
        vti_data.SetOrigin(-self._Lx/2, -self._Ly/2, -self._Lz/2)
        vti_data.SetSpacing(self._dx, self._dy, self._dz)

        n_data = vn.numpy_to_vtk(self._vals.reshape((self._Nx*self._Ny*self._Nz,3)))
        n_data.SetName("n")
        vti_data.GetPointData().AddArray(n_data)

        if self._mask_type is not None:
            mask_data = vn.numpy_to_vtk(
                self.mask_vals.reshape((self._Nx*self._Ny*self._Nz,1)))
            mask_data.SetName("lc_mask")
            vti_data.GetPointData().AddArray(mask_data)

        writer = vtkXMLImageDataWriter()
        writer.SetFileName(path)
        writer.SetInputData(vti_data)
        writer.Write()


    def get_pos(self, ix, iy, iz):
        """Returns the spatial position associated with the mesh indices (ix,iy,iz)

        It is assumed that the mesh is centered on the origin (0,0,0).
        """
        return (ix*self._dx-self._Lx/2, iy*self._dy-self._Ly/2, iz*self._dz-self._Lz/2)

    
    def get_mesh_dimensions(self):
        """Returns the dimensions (Nx,Ny,Nz) of the simulation mesh"""
        return (self._Nx, self._Ny, self._Nz)

    
    def get_mesh_lengths(self):
        """Returns the lengths (Lx,Ly,Lz) of the simulation mesh"""
        return (self._Lx, self._Ly, self._Lz)


    def get_mesh_spacings(self):
        """Returns the spacings (dx,dy,dz) of the simulation mesh"""
        return (self._dx, self._dy, self._dz)


    def get_n_vertices(self):
        """Returns the number of vertices in the simulation mesh"""
        return self._Nx*self._Ny*self._Nz
