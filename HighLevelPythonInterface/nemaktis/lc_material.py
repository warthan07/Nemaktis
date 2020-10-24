import numpy as np
import os

from vtk import vtkImageData, vtkXMLImageDataReader, vtkXMLImageDataWriter
from vtk.util import numpy_support as vn
from scipy.spatial.transform import Rotation as R
from scipy.interpolate import interpn

from warnings import simplefilter
simplefilter(action='ignore', category=FutureWarning)


class LCMaterial(object):
    """
    A class containing the LC orientational field data (director or q-tensor) and physics
    constants.

    Parameters
    ----------
    lc_field : :class:`~nemaktis.lc_material.DirectorField` or :class:`~nemaktis.lc_material.QTensorField` object
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
    def __init__(self, *, lc_field, ne, no, nhost = 1, nin = 1):
        self.lc_field = lc_field
        self.ne = ne
        self.no = no
        self.nhost = nhost
        self.nin = nin
        self.iso_layer_indices = [] 
        self.iso_layer_thicknesses = [] 


    def add_isotropic_layer(self, *, nlayer, thickness):
        """
        Add an isotropic layer above the sample.  Light is assumed to propagate in the
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



class TensorField:
    """The TensorField class stores the discrete data of a tensor field on a cartesian mesh.
    Currently, only first-order tensor (vector) fields and symmetric second-order tensor
    fields are supported by this class. The ordering of degree of freedoms for each mesh
    points is as follows:
    
    - Vector field n: (nx,ny,nz)
    - Symmetric second-order tensor field Q: (Qxx,Qyy,Qzz,Qxy,Qxz,Qyz)

    This class is initialized given either the lengths and dimensions of the associated 3D
    mesh or a path to a vti file containing the tensor field and mesh details.

    In the first version of this constructor:

    .. code-block:: python
    
        field = TensorField(
            mesh_lengths=(Lx,Ly,Lz), mesh_dimensions=(Nx,Ny,Nz), tensor_order=m)

    the actual values of the tensor field needs to be provided later via the setter method
    "vals", which should be given a numpy array of shape (Nz,Ny,Nx,Nv), where Nv=3 (Nv=6) if
    m=1 (m=2). The mesh lengths needs to be specified in micrometer.

    In the second version of this constructor:

    .. code-block:: python
    
        field = TensorField(vti_file="path to vti file", vti_array="name of tensor array")

    the values of the tensor field and the details of the mesh are automatically assigned
    from the given vti file and array name. 
    """

    _mask_type = None
    _mask_formula = None
    _mask_vals = None

    def __init__(self, **kwargs):
        if "vti_file" in kwargs and "vti_array" in kwargs:
            if not os.path.isfile(kwargs["vti_file"]):
                raise Exception("VTI file does not exists")

            print("{ Initializing tensor field from "+kwargs["vti_file"]+" }")

            reader = vtkXMLImageDataReader()
            reader.SetFileName(kwargs["vti_file"])
            reader.Update()

            vti_data = reader.GetOutput()
            dims = np.array(vti_data.GetDimensions())
            spacings = np.array(vti_data.GetSpacing())

            if not vti_data.GetPointData().HasArray(kwargs["vti_array"]):
                raise Exception( \
                    "Missing vtk array '%s' in the given vti file" % kwargs["vti_array"])

            tensor_data = \
                vn.vtk_to_numpy(vti_data.GetPointData().GetArray(kwargs["vti_array"]))
            field_dim = tensor_data.shape[1]
            if field_dim!=3 and field_dim!=6:
                raise Exception( \
                    "Unsupported tensor type, the field dimension should be 3 (vector field)"
                    "or 6 (symmetric second-order tensor field)")

            self._vals =  tensor_data.reshape((dims[2], dims[1], dims[0], field_dim))
            self._Nv = field_dim
            (self._Nx, self._Ny, self._Nz) = tuple(int(dim) for dim in dims)
            (self._dx, self._dy, self._dz) = tuple(spacings)
            (self._Lx, self._Ly, self._Lz) = tuple(spacings*(dims-1))

            if vti_data.GetPointData().HasArray("domain_mask"):
                mask_vals = vn.vtk_to_numpy(vti_data.GetPointData().GetArray("domain_mask"))
                self._mask_type = "raw"
                self._mask_vals = mask_vals.reshape((dims[2],dims[1],dims[0]))

        elif "mesh_lengths" in kwargs and "mesh_dimensions" in kwargs \
                and "tensor_order" in kwargs:
            if len(kwargs["mesh_lengths"])!=3:
                raise Exception("mesh_lengths should be an array-like object of size 3")
            if len(kwargs["mesh_dimensions"])!=3:
                raise Exception("mesh_dimensions should be an array-like object of size 3")

            print("{ Initializing empty tensor field }")

            dims = np.array(kwargs["mesh_dimensions"])
            lengths = np.array(kwargs["mesh_lengths"])

            if kwargs["tensor_order"]==1:
                self._Nv = 3
            elif kwargs["tensor_order"]==2:
                self._Nv = 6
            else:
                raise Exception(
                    "'tensor_order' should be 1 (vector field) or 2 (symmetric second-order "
                    "tensor field)")

            (self._Nx, self._Ny, self._Nz) = tuple(int(dim) for dim in dims)
            (self._Lx, self._Ly, self._Lz) = tuple(lengths)
            (self._dx, self._dy, self._dz) = tuple(lengths/np.maximum(np.ones(3),dims-1))
            self._vals = np.zeros((self._Nz,self._Ny,self._Nx,self._Nv))

        else:
            raise Exception("Could not parse the constructor parameters of TensorField")



    def set_mask(self, *,  mask_type, mask_formula = None, mask_ndarray = None):
        """Set a mask for the definition domain of this tensor field. This method allows to
        specifify an arbitrary complex definition domain which is a subset of the regular
        cartesian mesh specified at construction. Positive mask values are associated with
        the definition domain, while negative values are associated with the "host" embedding
        domain (for LC structure, this would be the host fluid or material which encases the
        LC domain). Three possible ways of initializing the mask are possible. If you simply
        want to specify a spherical domain for a droplet centered on the mesh and of
        diameter equal to the mesh length along z, call:

        .. code-block:: python
            
            tensor_field.set_mask(mask_type="droplet")

        You can also use a string formula depending on the space variables ``x``, ``y`` and
        ``z`` and which must evaluates to a value >=0 if the associated point is inside the
        definition domain of the tensor field, else to a value <=0:

        .. code-block:: python
            
            tensor_field.set_mask(mask_type="formula", mask_formula="your formula")

        Finally, you can directly gives a numpy array of shape (Nz,Ny,Nx), where each value
        in this array must be >=0 if the associated mesh point is inside the definition
        domain, else <=0:
        
        .. code-block:: python
            
            tensor_field.set_mask(mask_type="raw", mask_ndarray=your_mask_array)

        """
        if mask_type=="droplet":
            self._mask_type = mask_type
            self._mask_formula = str((self._Lz/2)**2)+"-x**2-y**2-z**2"
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
        """Returns the mask boolean array."""
        if self._mask_type=="raw":
            return self._mask_vals
        elif self._mask_type=="formula" or self._mask_type=="droplet":
            z,y,x = np.meshgrid(
                np.linspace(-0.5*self._dz*(self._Nz-1),0.5*self._dz*(self._Nz-1),self._Nz),
                np.linspace(-0.5*self._dy*(self._Ny-1),0.5*self._dy*(self._Ny-1),self._Ny),
                np.linspace(-0.5*self._dx*(self._Nx-1),0.5*self._dx*(self._Nx-1),self._Nx),
                indexing="ij")
            return eval(self._mask_formula)
        else:
            return None

    
    def extend(self, scale_x, scale_y):
        """
        Extend the computational mesh in the ``xy`` plane by padding new points near the ``x``
        and ``y`` boundaries. The associated new data points are initialized with the edge
        value of the tensor field on the ``x`` and ``y`` boundaries.

        Parameters
        ----------
        scale_x : float
            The mesh length in the x-direction will be scaled by this factor.
        scale_y : float
            The mesh length in the y-direction will be scaled by this factor.
        """
        print("{ Extending computational mesh }")
        
        if self._mask_vals is not None:
            raise Exception("You should always call 'set_mask' after 'extend', not before")

        Nx_pad = int(self._Nx*(scale_x-1)/2)
        Ny_pad = int(self._Ny*(scale_y-1)/2)
        self._vals = np.pad(
            self._vals, ((0,0),(Ny_pad,Ny_pad),(Nx_pad,Nx_pad),(0,0)), "edge")

        self._Nx = self._Nx+2*Nx_pad
        self._Ny = self._Ny+2*Ny_pad

        self._Lx = self._dx*(self._Nx-1)
        self._Ly = self._dy*(self._Ny-1)


    def rotate_90deg(self, axis):
        """
        Apply a solid rotation of the tensor field of 90 degrees around the specified axis.
        This is a lossless operation which does not rely on interpolation.
        
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

        # axis that will be flipped and components that will be reversed after rotation
        if len(axis)==1 or (len(axis)==2 and axis[0]=="+"):
            flip_ax = (ax+1)%3
            if self._Nv==3:
                negate_comps = [flip_ax]
            if self._Nv==6:
                negate_comps = [5-ax,3+(3-ax)%3]
        elif len(axis)==2 and axis[0]=="-":
            flip_ax = (ax+2)%3
            if self._Nv==3:
                negate_comps = [flip_ax]
            if self._Nv==6:
                negate_comps = [5-ax,3+(4-ax)%3]
        else:
            raise Exception("Could not parse axis.")

        print("{ 90 deg rotation around axis "+axis[-1]+" }")
        
        if self._mask_vals is not None:
            raise Exception("You should always call set_mask after rotate_90deg")

        # mapping for the axis order in the raw data array (axis 3 is field component axis)
        ax_map = [2, 1, 0, 3]

        # permutation of axes and components equivalent to the rotation operation
        permut_ax = np.arange(0,4)
        permut_ax[ax_map[(ax+1)%3]] = ax_map[(ax+2)%3]
        permut_ax[ax_map[(ax+2)%3]] = ax_map[(ax+1)%3]

        permut_comp = np.arange(0,self._Nv)
        permut_comp[(ax+1)%3] = (ax+2)%3
        permut_comp[(ax+2)%3] = (ax+1)%3
        if self._Nv==6:
            permut_comp[3+(3-ax)%3] = 3+(4-ax)%3
            permut_comp[3+(4-ax)%3] = 3+(3-ax)%3

        # We apply the axis and component permutations, followed by component negations
        self._vals = np.flip(
            np.transpose(self.vals, tuple(permut_ax)),
            axis = ax_map[flip_ax])[:,:,:,permut_comp]
        for comp in negate_comps:
            self._vals[:,:,:,comp] = -self._vals[:,:,:,comp]

        (self._Nx, self._Ny, self._Nz) = tuple(
            self.get_mesh_dimensions()[i] for i in permut_comp[0:3])
        (self._Lx, self._Ly, self._Lz) = tuple(
            self.get_mesh_lengths()[i] for i in permut_comp[0:3])
        (self._dx, self._dy, self._dz) = tuple(
            self.get_mesh_spacings()[i] for i in permut_comp[0:3])


    def rotate_180deg(self, axis):
        """
        Apply a solid rotation of the tensor field of 180 degrees around the specified axis.
        This is a lossless operation which does not rely on interpolation.
        
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
        
        if self._mask_vals is not None:
            raise Exception("You should always call set_mask after rotate_180deg")

        # mapping for the axis order in the raw data array (axis 3 is component axis)
        ax_map = np.array([2, 1, 0, 3])

        # Axes that will be flipped and components that will be reversed after rotation
        flip_axes = [(ax+1)%3, (ax+2)%3]
        if self._Nv==3:
            negate_comps = flip_axes
        if self._Nv==6:
            negate_comps = [3+(3-ax)%3, 3+(4-ax)%3]

        # We apply the rotation
        self._vals = np.flip(self.vals, axis = tuple(ax_map[flip_axes]))
        for comp in negate_comps:
            self._vals[:,:,:,comp] = -self._vals[:,:,:,comp]
        
    def rotate(self, axis, angle, fill_value=None):
        """
        Apply a solid rotation of the tensor field of an arbitrary angle around the
        specified axis. This is a lossy operation that will rely on interpolation, so
        possible artefacts can appear if the tensor field data is not smooth enough.
        
        Parameters
        ----------
        axis : str
            Axis around which to perform the rotation. Need to be under the form 'A' where
            'A'='x', 'y' or 'z' defines the rotation axis.
        angle : float
            Angle of rotation in degrees.
        """
        if axis[0]=="x":
            ax = 0
        elif axis[0]=="y":
            ax = 1
        elif axis[0]=="z":
            ax = 2
        else:
            raise Exception("Could not parse axis.")

        print("{ Rotation of %.2f° around axis %s }" % (angle,axis[0]))
        
        if self._mask_vals is not None:
            raise Exception("You should always call 'set_mask' after 'rotate', not before")
        
        u = np.zeros(3)
        u[ax] = 1
        rot_mat_inv = R.from_rotvec(-angle*np.pi/180*u).as_dcm()

        # For vector field, the transformation operator is simply the rotation matrix. For
        # symmetric second-order tensor field, the transformation operator can be obtained
        # in Mathematica with appropriate cartesian product and slicing operations (which we
        # compact based on roll and flip matrix operations)
        if self._Nv==3:
            transf_op = R.from_rotvec(angle*np.pi/180*u).as_dcm()
        if self._Nv==6:
            G = R.from_rotvec(angle*np.pi/180*u).as_dcm()
            transf_op = np.zeros((6,6))
            transf_op[0:3,0:3] = np.power(G,2)
            transf_op[0:3,3:6] = 2*np.flip(np.roll(G,1,axis=1)*np.roll(G,2,axis=1),axis=1)
            transf_op[3:6,0:3] = np.flip(np.roll(G,1,axis=0)*np.roll(G,2,axis=0),axis=0)
            transf_op[3:6,3:6] = np.flip(
                np.roll(np.roll(G,1,axis=1),2,axis=0)*np.roll(np.roll(G,2,axis=1),1,axis=0)+
                np.roll(np.roll(G,1,axis=1),1,axis=0)*np.roll(np.roll(G,2,axis=1),2,axis=0))
            
        x = np.linspace(-self._Lx/2, self._Lx/2, self._Nx)
        y = np.linspace(-self._Ly/2, self._Ly/2, self._Ny)
        z = np.linspace(-self._Lz/2, self._Lz/2, self._Nz)
        Z,Y,X = np.meshgrid(z,y,x,indexing="ij")
        
        pos = np.stack((X.flatten(),Y.flatten(),Z.flatten()), axis=1)
        pos_rot = np.dot(rot_mat_inv,pos.transpose()).transpose()
 
        tmp = interpn((z,y,x), self._vals, np.flip(pos_rot, axis=1),
                      bounds_error=False, fill_value=fill_value)
        self._vals = np.dot(transf_op,tmp.transpose()).transpose().reshape(
            (self._Nz,self._Ny,self._Nx,self._Nv))
        
        
    def rescale_mesh(self, scaling_factor):
        """
        Uniformly scale the mesh using the given scaling factor.
        
        Parameters
        ----------
        scaling_factor : factor
            The mesh lengths and spacings will be multiplied by this factor.
        """
        (self._Lx,self._Ly,self._Lz) = tuple(scaling_factor*np.array(self.get_mesh_lengths()))
        (self._dx,self._dy,self._dz) = tuple(scaling_factor*np.array(self.get_mesh_spacings()))
        

    @property
    def vals(self):
        """Numpy array for the tensor values, of shape (Nz,Ny,Nx,Nv), where Nv=3 for a vector
        field and Nv=6 for a symmetric second-order tensor field (we only store the
        [xx,yy,zz,xy,xz,yz] components for efficiency reasons).
        """
        return self._vals


    @vals.setter
    def vals(self, tensor_ndarray):
        if self._vals.shape==tensor_ndarray.shape:
            self._vals = tensor_ndarray
        else:
            raise Exception("Wrong shape for the tensor field ndarray")


    def save_to_vti(self, file_name, array_name):
        """Save the tensor field inside a vti file.

        Parameters
        ----------
        file_name : string 
            Path to the exported vti file. The ".vti" extension is automatically appended,
            no need to include it in this parameter (but in case you do only one extension
            will be added).
        array_name : string
            Name of the vti array that will store the tensor field.
        """ 
        if file_name[-4:]==".vti":
            path = file_name
        else:
            path = file_name+".vti"

        print("{ Saving tensor field to "+path+" }")
        vti_data = vtkImageData()
        vti_data.SetDimensions(self._Nx, self._Ny, self._Nz)
        vti_data.SetOrigin(-self._Lx/2, -self._Ly/2, -self._Lz/2)
        vti_data.SetSpacing(self._dx, self._dy, self._dz)

        tensor_data = \
            vn.numpy_to_vtk(self._vals.reshape((self._Nx*self._Ny*self._Nz,self._Nv)))
        tensor_data.SetName(array_name)
        vti_data.GetPointData().AddArray(tensor_data)

        if self._mask_type is not None:
            mask_data = vn.numpy_to_vtk(
                self.mask_vals.reshape((self._Nx*self._Ny*self._Nz,1)))
            mask_data.SetName("domain_mask")
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


class DirectorField(TensorField):
    """
    A specialization of the TensorField class for director fields.

    The two versions of the constructor of the parent class TensorField are simplified since
    we do not need the parameters 'tensor_order' (always 1 for a director field) or
    'vti_array' (assumed to be "n"):

    .. code-block:: python
    
        # First version of the constructor
        nfield = DirectorField(
            mesh_lengths=(Lx,Ly,Lz), mesh_dimensions=(Nx,Ny,Nz))
        # Second version of the constructor
        nfield = DirectorField(
            vti_file="path to vti file", vti_array="name of tensor array")

    In addition to all the methods of the parent class for initializing and manipulating 
    the field values, we specialize the "save_to_vti" method (imposing that the exported vti
    array name is always "n") and provide additional methods for initializing the
    director field values from theoretical functions, exporting a q-tensor field from the
    director field, and normalizing the director field to unit norm.
    """
    def __init__(self, **kwargs):
        if "vti_file" in kwargs:
            kwargs["vti_array"] = "n"
        elif "mesh_lengths" in kwargs and "mesh_dimensions" in kwargs:
            kwargs["tensor_order"] = 1
        else:
            raise Exception("Could not parse the constructor parameters of DirectorField")
        super().__init__(**kwargs)
    

    def init_from_funcs(self, nx_func, ny_func, nz_func):
        """Initialize the director field from three functions for each of its component. The
        functions must depend on the space variables ``x``, ``y`` and ``z``. We recall that
        the mesh is centered on the origin.

        If the given functions are numpy-vectorizable, this function should be pretty fast. If
        not, a warning will be printed and the faulty function(s) will be vectorized with the
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
        """Normalize the director field values to unit norm."""

        print("{ Normalizing director values to 1 }")
        norms = np.sqrt(np.sum(self._vals**2, axis=3, keepdims=True))
        norms[norms==0] = 1
        self._vals = self._vals / np.tile(norms, (1,1,1,3))


    def get_qtensor_field(self):
        """Returns a QTensorField object equivalent to the director field represented by this
        class, assuming a constant scalar order parameter (equal to its equilibrium value).

        Since in Nemaktis a q-tensor field is always renormalized by the equilibrium value
        of the order parameter, this method simply uses the formula Q_ij=(3*n_i*n_j-d_ij)/2
        to convert a director value n to a q-tensor value Q, with d_ij the kronecker delta. 
        """
        qfield = QTensorField(
            mesh_dimensions=self.get_mesh_dimensions(),
            mesh_lengths=self.get_mesh_lengths())

        qfield._vals[:,:,:,0] = (3*self._vals[:,:,:,0]**2-1)/2
        qfield._vals[:,:,:,1] = (3*self._vals[:,:,:,1]**2-1)/2
        qfield._vals[:,:,:,2] = (3*self._vals[:,:,:,2]**2-1)/2
        qfield._vals[:,:,:,3] = 3*self._vals[:,:,:,0]*self._vals[:,:,:,1]/2
        qfield._vals[:,:,:,4] = 3*self._vals[:,:,:,0]*self._vals[:,:,:,2]/2
        qfield._vals[:,:,:,5] = 3*self._vals[:,:,:,1]*self._vals[:,:,:,2]/2

        return qfield


    def save_to_vti(self, file_name):
        """Save the director field into a vti file, assuming "n" for the vti array name.

        Parameters
        ----------
        file_name : string 
            Path to the exported vti file. The ".vti" extension is automatically appended,
            no need to include it in this parameter (but in case you do only one extension
            will be added).
        """ 
        super().save_to_vti(file_name,"n")


class QTensorField(TensorField):
    """
    A specialization of the TensorField class for Q-tensor fields (i.e. the full LC order
    parameter field). 

    In Nemaktis, we always assume that the Q-tensor field is renormalized by the equilibrium
    value S_eq of the scalar order parameter S. This means for uniaxial LC, the Q-tensor
    far from topological defects can always be put under the following form:

    .. math::

        Q_{ij} = \\left(3n_in_j-\delta_{ij}\\right)/2

    The two versions of the constructor of the parent class TensorField are simplified since
    we do not need the parameters 'tensor_order' (always 2 for Q-tensor) or 'vti_array'
    (assumed to be "Q"):

    .. code-block:: python
    
        # First version of the constructor
        qfield = QTensorField(
            mesh_lengths=(Lx,Ly,Lz), mesh_dimensions=(Nx,Ny,Nz))
        # Second version of the constructor
        qfield = QTensorField(
            vti_file="path to vti file", vti_array="name of tensor array")

    In addition to all the methods of the parent class for initializing and manipulating 
    the field values, we specialize the "save_to_vti" method (imposing that the exported vti
    array name is always "Q") and provide addional methods for initializing the
    Q-tensor field values from theoretical functions, exporting a director field from a
    q-tensor field, and imposing the traceless constraint Tr(Q)=0.
    """
    def __init__(self, **kwargs):
        if "vti_file" in kwargs:
            kwargs["vti_array"] = "Q"
        elif "mesh_lengths" in kwargs and "mesh_dimensions" in kwargs:
            kwargs["tensor_order"] = 2
        else:
            raise Exception("Could not parse the constructor parameters of QTensorField")
        super().__init__(**kwargs)


    def init_from_funcs(self, Qxx_func, Qyy_func, Qzz_func, Qxy_func, Qxz_func, Qyz_func):
        """Initialize the Q-tensor field from six functions for each of its component (xx,
        yy, zz, xy, xz, yz). The functions must depend on the space variables ``x``, ``y``
        and ``z``. We recall that the mesh is centered on the origin.

        If the given functions are numpy-vectorizable, this function should be pretty fast.
        If not, a warning will be printed and the faulty function(s) will be vectorized with the
        numpy method ``vectorize`` (in which case you should expect a much slower execution
        time).
        """

        print("{ Calculating Q-tensor values from user functions }")

        # We verify if the user functions are vectorizable
        dummy_arr = np.ones((2,2,2))
        try:
            Qxx_func(dummy_arr,dummy_arr,dummy_arr)
        except:
            print("\tQxx_func is not vectorized, using a non-optimized version instead.")
            Qxx_func = np.vectorize(nz_func)
        try:
            Qyy_func(dummy_arr,dummy_arr,dummy_arr)
        except:
            print("\tQyy_func is not vectorized, using a non-optimized version instead.")
            Qyy_func = np.vectorize(nz_func)
        try:
            Qzz_func(dummy_arr,dummy_arr,dummy_arr)
        except:
            print("\tQzz_func is not vectorized, using a non-optimized version instead.")
            Qzz_func = np.vectorize(nz_func)
        try:
            Qxy_func(dummy_arr,dummy_arr,dummy_arr)
        except:
            print("\tQxy_func is not vectorized, using a non-optimized version instead.")
            Qxy_func = np.vectorize(nz_func)
        try:
            Qxz_func(dummy_arr,dummy_arr,dummy_arr)
        except:
            print("\tQxz_func is not vectorized, using a non-optimized version instead.")
            Qxz_func = np.vectorize(nz_func)
        try:
            Qyz_func(dummy_arr,dummy_arr,dummy_arr)
        except:
            print("\tQyz_func is not vectorized, using a non-optimized version instead.")
            Qyz_func = np.vectorize(nz_func)

        zz, yy, xx = np.meshgrid(np.linspace(-self._Lz/2, self._Lz/2, self._Nz),
                                 np.linspace(-self._Ly/2, self._Ly/2, self._Ny),
                                 np.linspace(-self._Lx/2, self._Lx/2, self._Nx),
                                 indexing="ij")
        self._vals = np.concatenate((np.expand_dims(Qxx_func(xx, yy, zz), axis=3),
                                     np.expand_dims(Qyy_func(xx, yy, zz), axis=3),
                                     np.expand_dims(Qzz_func(xx, yy, zz), axis=3),
                                     np.expand_dims(Qxy_func(xx, yy, zz), axis=3),
                                     np.expand_dims(Qxz_func(xx, yy, zz), axis=3),
                                     np.expand_dims(Qyz_func(xx, yy, zz), axis=3)), 3)


    def apply_traceless_constraint(self):
        """Apply the traceless constraint Tr(Q)=0 to this q-tensor field."""
        self._vals[:,:,:,0:3] -= np.sum(self._vals[:,:,:,0:3],axis=-1,keepdims=True)/3


    def get_director_field(self):
        """Returns a DirectorField object equivalent to the q-tensor field represented by this
        class, assuming a unixial medium and discarding any variations of the scalar order
        parameter S.

        In practice this method simply calculate the eigenvectors of the Q-tensor field and
        initialize the director field from the eigenvectors with highest algebraic value.
        The returned director field is therefore not fully equivalent to the Q-tensor field
        if there are topological defects or biaxiality inside the LC structure.
        """

        # First, we calculate analytically the eigenvalues
        # (algorithm from wikipedia page "Eigenvalue algorithm",
        # simplified a bit because Tr[Q]=0)
        qxx = self._vals[:,:,:,0]
        qyy = self._vals[:,:,:,1]
        qzz = self._vals[:,:,:,2]
        qxy = self._vals[:,:,:,3]
        qxz = self._vals[:,:,:,4]
        qyz = self._vals[:,:,:,5]

        p1 = qxy**2 + qxz**2 + qyz**2
        p2 = np.sqrt((qxx**2 + qyy**2 + qzz**2 + 2*p1)/6.)
        det = (qxx*qyy*qzz + 2*qxy*qyz*qxz - qxz*qyy*qxz - qxy*qxy*qzz - qyz*qyz*qxx) / p2**3
        phi = np.arccos(np.clip(det/2,-1,1))/3

        eig1 = 2*p2*np.cos(phi)
        eig3 = 2*p2*np.cos(phi+2*np.pi/3)
        eig2 = -eig1-eig3


        # Then, we calculate the director as the unit eigenvector associated with eig1.
        # Since Q is symmetric (and therefore normal), we can calculate this eigenvector
        # with a cross product between any two non-colinear column vectors of Q-eig1*I.
        # I do not want to waste time writing a proper numba vectorized algorithm for
        # this, so I simply calculate all possible candidates for the eigenvectors and
        # select the ones with highest (non-zero) norms, which can simply done using
        # native broadcasting rules of numpy.
        eigenvecs = np.zeros(qxx.shape+(3,3))
        eigenvecs[:,:,:,0,0] = qxy*qyz - qxz*(qyy-eig1)
        eigenvecs[:,:,:,0,1] = qxy*qxz - qyz*(qxx-eig1)
        eigenvecs[:,:,:,0,2] = (qxx-eig1)*(qyy-eig1) - qxy*qxy

        eigenvecs[:,:,:,1,0] = (qyy-eig1)*(qzz-eig1) - qyz*qyz
        eigenvecs[:,:,:,1,1] = qxz*qyz - (qzz-eig1)*qxy
        eigenvecs[:,:,:,1,2] = qxy*qyz - (qyy-eig1)*qxz

        eigenvecs[:,:,:,2,0] = qxy*(qzz-eig1) - qxz*qyz
        eigenvecs[:,:,:,2,1] = qxz*qxz - (qzz-eig1)*(qxx-eig1)
        eigenvecs[:,:,:,2,2] = (qxx-eig1)*qyz - qxy*qxz

        norms = np.sqrt(np.sum(eigenvecs**2,axis=-1))
        max_ids = np.argmax(norms, axis=-1)[:,:,:,np.newaxis,np.newaxis]

        nfield = DirectorField(
            mesh_dimensions=self.get_mesh_dimensions(),
            mesh_lengths=self.get_mesh_lengths())
        nfield.vals = np.squeeze( \
            np.take_along_axis(eigenvecs, max_ids, axis=3)
            / np.take_along_axis(norms[:,:,:,:,np.newaxis], max_ids, axis=3))

        return nfield


    def save_to_vti(self, file_name):
        """Save the Q-tensor field into a vti file, assuming "Q" for the vti array name.

        Parameters
        ----------
        file_name : string 
            Path to the exported vti file. The ".vti" extension is automatically appended,
            no need to include it in this parameter (but in case you do only one extension
            will be added).
        """ 
        super().save_to_vti(file_name,"Q")
