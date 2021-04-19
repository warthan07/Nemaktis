import numpy as np
import os

from vtk import vtkImageData, vtkXMLImageDataReader, vtkXMLImageDataWriter, vtkStringArray
from vtk.util import numpy_support as vn
from scipy.spatial.transform import Rotation as R
from scipy.interpolate import interpn

from warnings import simplefilter
simplefilter(action='ignore', category=FutureWarning)


class TensorField:
    """The TensorField class stores the discrete data of a tensor field on a cartesian mesh.
    Currently, only zeroth-order (scalar) fields, first-order tensor (vector) fields and
    symmetric second-order tensor fields are supported by this class. For the last two
    cases, the ordering of degree of freedoms for each mesh points is as follows:
    
    - Vector field n: (nx,ny,nz)
    - Symmetric second-order tensor field Q: (Qxx,Qyy,Qzz,Qxy,Qxz,Qyz)

    This class is initialized given either a path to a vti file containing the tensor field
    and mesh details, or the lengths and dimensions of the associated 3D mesh, the tensor
    order and the name of the tensor (for logging and export purpose)

    In the first version of this constructor:

    .. code-block:: python
    
        nfield = TensorField(vti_file="path to vti file", vti_array="name of tensor array")

    the values of the tensor field and the details of the mesh are automatically assigned
    from the given vti file and array name (the latter also being used as tensor name).

    In the second version of this constructor:

    .. code-block:: python
    
        nfield = TensorField(
            mesh_lengths=(Lx,Ly,Lz), mesh_dimensions=(Nx,Ny,Nz),
            tensor_order=m, tensor_name=name)

    the actual values of the tensor field needs to be provided later via the setter method
    "vals", which should be given a numpy array of shape (Nz,Ny,Nx,Nv), where Nv=1,3,6 if
    m=0,1,2. The mesh lengths needs to be specified in micrometer.
    """
    def __init__(self, **kwargs):
        if "vti_file" in kwargs and "vti_array" in kwargs:
            if not os.path.isfile(kwargs["vti_file"]):
                raise Exception("VTI file does not exists")

            self._name = kwargs["vti_array"]
            print("{ Initializing tensor field %s from %s }" \
                % (self._name,kwargs["vti_file"]))

            reader = vtkXMLImageDataReader()
            reader.SetFileName(kwargs["vti_file"])
            reader.Update()

            vti_data = reader.GetOutput()
            dims = np.array(vti_data.GetDimensions())
            spacings = np.array(vti_data.GetSpacing())

            if not vti_data.GetPointData().HasArray(self._name):
                raise Exception( \
                    "Missing vtk array '%s' in the given vti file" % self._name)

            tensor_data = \
                vn.vtk_to_numpy(vti_data.GetPointData().GetArray(self._name))
            field_dim = \
                vti_data.GetPointData().GetArray(self._name).GetNumberOfComponents()
            if field_dim!=1 and field_dim!=3 and field_dim!=6:
                raise Exception( \
                    "Unsupported tensor type, the field dimension should be 1 (scalar "
                    "field), 3 (vector field) or 6 (symmetric second-order tensor field)")

            self._vals =  tensor_data.reshape((dims[2], dims[1], dims[0], field_dim))
            self._Nv = field_dim
            (self._Nx, self._Ny, self._Nz) = tuple(int(dim) for dim in dims)
            (self._dx, self._dy, self._dz) = tuple(spacings)
            (self._Lx, self._Ly, self._Lz) = tuple(spacings*(dims-1))

        elif "mesh_lengths" in kwargs and "mesh_dimensions" in kwargs \
                and "tensor_order" in kwargs and "tensor_name" in kwargs:
            if len(kwargs["mesh_lengths"])!=3:
                raise Exception("mesh_lengths should be an array-like object of size 3")
            if len(kwargs["mesh_dimensions"])!=3:
                raise Exception("mesh_dimensions should be an array-like object of size 3")

            self._name = kwargs["tensor_name"]
            print("{ Initializing tensor field %s with zero values }" % self._name)

            dims = np.array(kwargs["mesh_dimensions"])
            lengths = np.array(kwargs["mesh_lengths"])

            if kwargs["tensor_order"]==0:
                self._Nv = 1
            elif kwargs["tensor_order"]==1:
                self._Nv = 3
            elif kwargs["tensor_order"]==2:
                self._Nv = 6
            else:
                raise Exception(
                    "'tensor_order' should be 0 (scalar field), 1 (vector field) or 2 "
                    "(symmetric second-order tensor field)")

            (self._Nx, self._Ny, self._Nz) = tuple(int(dim) for dim in dims)
            (self._Lx, self._Ly, self._Lz) = tuple(lengths)
            (self._dx, self._dy, self._dz) = tuple(lengths/np.maximum(np.ones(3),dims-1))
            self._vals = np.zeros((self._Nz,self._Ny,self._Nx,self._Nv))

        else:
            raise Exception("Could not parse the constructor parameters of TensorField")

    
    def extend(self, scale_x, scale_y):
        """
        Extend the computational mesh in the ``xy`` plane by padding new points near the
        ``x`` and ``y`` boundaries. The associated new data points are initialized with the
        edge value of the tensor field on the ``x`` and ``y`` boundaries.

        Parameters
        ----------
        scale_x : float
            The mesh length in the x-direction will be scaled by this factor.
        scale_y : float
            The mesh length in the y-direction will be scaled by this factor.
        """
        print("{ Extending computational mesh of tensor field %s }" % self._name)

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
            Axis around which to perform the rotation. Need to be under the form '[s]A'
            where the optional parameter 's'='+' or '-' decribes the sign of rotation and
            'A'='x', 'y' or 'z' defines the rotation axis.
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
            if self._Nv==1:
                negate_comps = []
            if self._Nv==3:
                negate_comps = [flip_ax]
            if self._Nv==6:
                negate_comps = [5-ax,3+(3-ax)%3]
        elif len(axis)==2 and axis[0]=="-":
            flip_ax = (ax+2)%3
            if self._Nv==1:
                negate_comps = []
            if self._Nv==3:
                negate_comps = [flip_ax]
            if self._Nv==6:
                negate_comps = [5-ax,3+(4-ax)%3]
        else:
            raise Exception("Could not parse axis.")

        print("{ 90 deg rotation of tensor field %s around axis %s }" \
            % (self._name,axis[-1]))

        # mapping for the axis order in the raw data array (axis 3 is field component axis)
        ax_map = [2, 1, 0, 3]

        # permutation of axes and components equivalent to the rotation operation
        permut_ax = np.arange(0,4)
        permut_ax[ax_map[(ax+1)%3]] = ax_map[(ax+2)%3]
        permut_ax[ax_map[(ax+2)%3]] = ax_map[(ax+1)%3]

        permut_comp = np.arange(0,self._Nv)
        if self._Nv!=1:
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
            self.get_mesh_dimensions()[i] for i in permut_ax[0:3])
        (self._Lx, self._Ly, self._Lz) = tuple(
            self.get_mesh_lengths()[i] for i in permut_ax[0:3])
        (self._dx, self._dy, self._dz) = tuple(
            self.get_mesh_spacings()[i] for i in permut_ax[0:3])


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

        print("{ 180 deg rotation of tensor field %s around axis %s }" \
            % (self._name,axis[-1]))

        # mapping for the axis order in the raw data array (axis 3 is component axis)
        ax_map = np.array([2, 1, 0, 3])

        # Axes that will be flipped and components that will be reversed after rotation
        flip_axes = [(ax+1)%3, (ax+2)%3]
        if self._Nv==1:
            negate_comps = []
        if self._Nv==3:
            negate_comps = flip_axes
        if self._Nv==6:
            negate_comps = [3+(3-ax)%3, 3+(4-ax)%3]

        # We apply the rotation
        self._vals = np.flip(self.vals, axis = tuple(ax_map[flip_axes]))
        for comp in negate_comps:
            self._vals[:,:,:,comp] = -self._vals[:,:,:,comp]

        
    def rotate(self, axis, angle, *, method="linear", fill_value=None):
        """
        Apply a solid rotation of the tensor field of an arbitrary angle around the
        specified axis. This is a lossy operation that will rely on interpolation, so
        possible artefacts can appear if the tensor field data is not smooth enough.
        
        Note that the computational box stays unchanged after calling this method, which
        means that all mesh points falling out of the computational box after rotation will
        be clipped.  Please be aware that by default, new tensor values near the mesh
        boundaries will possibly be interpolated from extrapolated values outside the mesh
        since a regular cartesian mesh is not rotationally invariant. These values are
        generally spurious, and should always be ignored through the use of a domain mask
        fully included inside the cartesian mesh whatever the applied rotation (e.g. a
        sphere centered inside the mesh). Alternatively, you can specify the fill_value
        parameter, which should contain the (uniform) value of the tensor field for all
        points outside the cartesian mesh.
        
        Parameters
        ----------
        axis : str
            Axis around which to perform the rotation. Need to be under the form 'A' where
            'A'='x', 'y' or 'z' defines the rotation axis.
        angle : float
            Angle of rotation in degrees.
        method : "linear" (default) or "nearest", optional named parameter
            Interpolation method when rotating the mesh.
        fill_value : array, optional named parameter
            Value of the tensor field for all points outside the cartesian mesh.
        """
        if axis[0]=="x":
            ax = 0
        elif axis[0]=="y":
            ax = 1
        elif axis[0]=="z":
            ax = 2
        else:
            raise Exception("Could not parse axis.")

        print("{ %.2f deg rotation of tensor field %s around axis %s }" \
            % (angle,self._name,axis[0]))
        
        u = np.zeros(3)
        u[ax] = 1
        rot_mat_inv = R.from_rotvec(-angle*np.pi/180*u).as_dcm()

        # For vector field, the transformation operator is simply the rotation matrix. For
        # symmetric second-order tensor field, the transformation operator can be obtained
        # in Mathematica with appropriate cartesian product and slicing operations (which we
        # compact based on roll and flip matrix operations). For scalar field, no need to
        # apply a transformation operator since only the mesh is modified in this case.
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
        if self._Nv==1:
            self._vals = tmp.reshape((self._Nz,self._Ny,self._Nx,self._Nv))
        else:
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

        print("{ Rescaling mesh of tensor field %s by a factor %.1f }" \
            % (self._name,scaling_factor))

        (self._Lx,self._Ly,self._Lz) = tuple(
            scaling_factor*np.array(self.get_mesh_lengths()))
        (self._dx,self._dy,self._dz) = tuple(
            scaling_factor*np.array(self.get_mesh_spacings()))
        

    @property
    def vals(self):
        """Numpy array for the tensor values, of shape (Nz,Ny,Nx,Nv), where Nv=1 for a
        scalar field, Nv=3 for a vector field and Nv=6 for a symmetric second-order tensor
        field (we only store the [xx,yy,zz,xy,xz,yz] components for efficiency reasons).
        """
        return self._vals


    @vals.setter
    def vals(self, tensor_ndarray):
        if self._vals.shape==tensor_ndarray.shape:
            self._vals = tensor_ndarray
        else:
            raise Exception("Wrong shape for the tensor field ndarray")


    def save_to_vti(self, file_name):
        """Save the tensor field inside a vti file. The name of the tensor field given at
        construction will be used as vti array name.

        Parameters
        ----------
        file_name : string 
            Path to the exported vti file. The ".vti" extension is automatically appended,
            no need to include it in this parameter (but in case you do only one extension
            will be added).
        """ 
        if file_name[-4:]==".vti":
            path = file_name
        else:
            path = file_name+".vti"

        print("{ Saving tensor field %s to %s }" % (self._name,path))
        vti_data = vtkImageData()
        vti_data.SetDimensions(self._Nx, self._Ny, self._Nz)
        vti_data.SetOrigin(-self._Lx/2, -self._Ly/2, -self._Lz/2)
        vti_data.SetSpacing(self._dx, self._dy, self._dz)

        tensor_data = \
            vn.numpy_to_vtk(self._vals.reshape((self._Nx*self._Ny*self._Nz,self._Nv)))
        tensor_data.SetName(self._name)
        vti_data.GetPointData().AddArray(tensor_data)

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
    A specialization of the TensorField class for director fields (i.e. unit vector fields)

    The constructors of the parent class TensorField are simplified since we do not need the
    parameter 'tensor_order' (always 1 for a director field):

    .. code-block:: python
    
        # First version of the constructor
        nfield = DirectorField(
            vti_file="path to vti file", vti_array="name of tensor array")
        # Second version of the constructor
        nfield = DirectorField(
            mesh_lengths=(Lx,Ly,Lz), mesh_dimensions=(Nx,Ny,Nz),
            tensor_name="name of tensor array")

    In addition to all the methods of the parent class for initializing and manipulating 
    the field values, we provide additional methods for initializing the
    director field values from theoretical functions, exporting a q-tensor field from the
    director field, and normalizing the director field to unit norm.
    """
    def __init__(self, **kwargs):
        if "mesh_lengths" in kwargs and "mesh_dimensions" in kwargs \
                and "tensor_name" in kwargs:
            kwargs["tensor_order"] = 1
        elif "vti_file" not in kwargs or "vti_array" not in kwargs:
            raise Exception("Could not parse the constructor parameters of DirectorField")
        super().__init__(**kwargs)

        # We check that this object is indeed a vector field in case of import from vti file
        if self._Nv!=3:
            raise Exception(
                "The given vti_array name does not correspond to a vector field")
    

    def init_from_funcs(self, nx_func, ny_func, nz_func):
        """Initialize the director field from three functions for each of its component. The
        functions must depend on the space variables ``x``, ``y`` and ``z``. We recall that
        the mesh is centered on the origin.

        If the given functions are numpy-vectorizable, this function should be pretty fast.
        If not, a warning will be printed and the faulty function(s) will be vectorized with
        the numpy method ``vectorize`` (in which case you should expect a much slower
        execution time).
        """

        print("{ Calculating values of director field %s from user functions }" \
            % self._name)
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

        print("{ Normalizing values of director field %s to unit norm }" % self._name)
        norms = np.sqrt(np.sum(self._vals**2, axis=3, keepdims=True))
        norms[norms==0] = 1
        self._vals = self._vals / np.tile(norms, (1,1,1,3))


    def get_qtensor_field(self):
        """Returns a QTensorField object equivalent to the director field represented by
        this class, assuming a constant scalar order parameter equal to its equilibrium
        value.

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
    we do not need the parameters 'tensor_order' (always 2 for Q-tensor) or
    'vti_array'/'tensor_name' (always assumed to be "Q"):

    .. code-block:: python
    
        # First version of the constructor
        nfield = QTensorField(
            vti_file="path to vti file", vti_array="name of tensor array")
        # Second version of the constructor
        nfield = QTensorField(
            mesh_lengths=(Lx,Ly,Lz), mesh_dimensions=(Nx,Ny,Nz))

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
            kwargs["tensor_name"] = "Q"
        else:
            raise Exception("Could not parse the constructor parameters of QTensorField")
        super().__init__(**kwargs)

        # We check that this object is indeed a vector field in case of import from vti file
        if self._Nv!=6:
            raise Exception(
                "The given vti_array name does not correspond to a symmetric second-order "
                "tensor")


    def init_from_funcs(self, Qxx_func, Qyy_func, Qzz_func, Qxy_func, Qxz_func, Qyz_func):
        """Initialize the Q-tensor field from six functions for each of its component (xx,
        yy, zz, xy, xz, yz). The functions must depend on the space variables ``x``, ``y``
        and ``z``. We recall that the mesh is centered on the origin.

        If the given functions are numpy-vectorizable, this function should be pretty fast.
        If not, a warning will be printed and the faulty function(s) will be vectorized with
        the numpy method ``vectorize`` (in which case you should expect a much slower
        execution time).
        """

        print("{ Calculating values of tensor field %s from user functions }" % self._name)

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

        print("{ Enforcing traceless constraint on tensor field %s }" % self._name)
        self._vals[:,:,:,0:3] -= np.sum(self._vals[:,:,:,0:3],axis=-1,keepdims=True)/3


    def get_director_field(self):
        """Returns a DirectorField object equivalent to the q-tensor field represented by
        this class, assuming a unixial medium and discarding any variations of the scalar
        order parameter S.

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
        det = (qxx*qyy*qzz + 2*qxy*qyz*qxz - qxz*qyy*qxz - qxy*qxy*qzz - qyz*qyz*qxx) \
              / p2**3
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


class MicroscopeSample(TensorField):
    """
    A class which stores all the information necessary to assemble the permittivity tensor
    field of a microscope sample. It allows to specify arbitrary number of subdomains on a
    regular cartesian mesh, all of which can be optically isotropic, unixial or biaxial, as
    well as any number of isotropic layers above or below the cartesian mesh (usually
    corresponding to sample plates confining the sample). The subdomain ids (default
    domain=0, additional subdomain ids=1,2,...) are internally stored in a scalar field
    defined on the same regular cartesian mesh representing the sample.
    
    Depending on the symmetry of the subdomains, you will need to specify a certain number
    of optical axis fields and refractive indices allowing to fully reconstruct the
    permittivity tensor during the microscope simulation. In the first version of the
    constructor, you need to specify the lengths and dimensions of the 3D regular cartesian
    mesh and an additional parameter allowing to specify the supported type of structures:

    .. code-block:: python
    
        # In this version of the constructor, no optical axis field need to be specified and
        # therefore only isotropic subdomains are supported.
        sample = MicroscopeSample(
            mesh_lengths=(Lx,Ly,Lz), mesh_dimensions=(Nx,Ny,Nz),
            optical_axis_field_type=None)

        # In this version of the constructor, one optical axis field (called the director
        # field "n") need to be specified and therefore only isotropic or unixial subdomains
        # are supported.
        sample = MicroscopeSample(
            mesh_lengths=(Lx,Ly,Lz), mesh_dimensions=(Nx,Ny,Nz),
            optical_axis_field_type="Director")

        # In this version of the constructor, two optical axis fields (called the director
        # fields "n" and "m") need to be specified and therefore all three types of optical
        # symmetry (isotropic, unixial and biaxial) are supported for the subdomains.
        sample = MicroscopeSample(
            mesh_lengths=(Lx,Ly,Lz), mesh_dimensions=(Nx,Ny,Nz),
            optical_axis_field_type="DualDirector")

        # In this version of the constructor, no optical axis field are specified and
        # instead we use the full tensorial order parameter "Q" of a liquid crystal phase to
        # directly initialize the permittivity tensor field (knowing extraordinary and
        # ordinary indices)
        sample = MicroscopeSample(
            mesh_lengths=(Lx,Ly,Lz), mesh_dimensions=(Nx,Ny,Nz),
            optical_axis_field_type="QTensor")

    After initializing the MicroscopeSample with one of the constructor above, you need to
    appropriately initialize the values of required optical axis fields (nfield in
    "Director" mode, nfield+mfield in "DualDirector" mode, qfield in "QTensorField"), which
    are stored as attributes of this class. Please refer to the
    :class:`~nemaktis.lc_material.DirectorField` and
    :class:`~nemaktis.lc_material.QTensorField` classes on how to do that. After this step,
    you need to specify any relevant subdomains (by default, the cartesian mesh is
    initialized with the domain id "0", so if you only have one domain you do not need to do
    anything) and associated refractive indices. 

    Alternatively, you can get a fully initialized MicroscopeSample by specifying a path to
    a vti file exported by a previous instance of this class:

    .. code-block:: python
    
        sample = MicroscopeSample(vti_file="path to vti file")

    This will import the subdomain id field and any relevant optical axis fields, all the
    associated refractive indices, and all the details of isotropic layers above and below
    the sample.
    """
    def __init__(self, **kwargs):
        if "vti_file" in kwargs:
            ########
            # TODO #
            ########
            print("not yet implemented")

        elif "mesh_lengths" in kwargs and "mesh_dimensions" in kwargs \
                and "optical_axis_field_type" in kwargs:
            kwargs["tensor_order"] = 0
            kwargs["tensor_name"] = "DomainID"
            super().__init__(**kwargs)

            if kwargs["optical_axis_field_type"] \
                    in [None,"Director","DualDirector","QTensor"]:
                self._optical_axis_field_type = kwargs["optical_axis_field_type"]
            else:
                raise Exception( \
                    "Unrecognised value for optical_axis_field_type")

            # We initialize relevant optical axis field(s)
            if self._optical_axis_field_type=="Director":
                kwargs["tensor_name"] = "n"
                self._nfield = DirectorField(**kwargs)
            elif self._optical_axis_field_type=="DualDirector":
                kwargs["tensor_name"] = "n"
                self._nfield = DirectorField(**kwargs)
                kwargs["tensor_name"] = "m"
                self._mfield = DirectorField(**kwargs)
            elif self._optical_axis_field_type=="QTensor":
                kwargs["tensor_name"] = "Q"
                self._qfield = QTensorField(**kwargs)

            # Default refractive indices are invalid ("") and need to be initialized
            # later
            self._refractive_indices = np.empty((1,3), dtype=object)
            self._refractive_indices[:] = ""

            # By default there are no isotropic layers
            self._up_iso_layer_indices = [] 
            self._up_iso_layer_thicknesses = [] 
            self._lo_iso_layer_indices = [] 
            self._lo_iso_layer_thicknesses = [] 

            # We impose that the subdomain id field is integer
            self._vals = self._vals.astype(np.int32)

        else:
            raise Exception(
                "Could not parse the constructor parameters of MicroscopeSample")


    def add_subdomain(self, *,  mask_type, mask_formula = None, mask_ndarray = None):
        """Specify a new subdomain inside the regular cartesian mesh of this object. This
        new subdomain will be attributed a new identifier N+1, where N is the current
        number of subdomains. The identifier 0 is reserved for the default domain
        encompassed by the regular cartesian mesh, and any subsequent subdomain is
        substracted from this default domain. 

        The new subdomain is specified through a mask field defined on the same cartesian
        mesh as this class. True values are associated with the new subdomain, while False
        values correspond to the default domain or other subdomain already specified. Note
        that we check that the new subdomain is not intersecting with previously specified
        subdomains (but of course it should intersect with the default domain!).

        Three possible ways of initializing the subdomain mask are possible. If you simply
        want to specify a spherical domain for a droplet centered on the mesh and of
        diameter equal to the mesh length along z, call:

        .. code-block:: python

            sample.add_subdomain(mask_type="droplet")

        You can also use a string formula depending on the space variables ``x``, ``y`` and
        ``z`` and which must evaluates to a value >=0 if the associated point is inside the
        new subdomain, else to a value <=0:

        .. code-block:: python

            sample.add_subdomain(mask_type="formula", mask_formula="your formula")
            
        Finally, you can directly gives a numpy array of shape (Nz,Ny,Nx), where each value
        in this array must be True if the associated mesh point is inside the new subdomain,
        else False:

        .. code-block:: python

            sample.add_subdomain(mask_type="raw", mask_ndarray=your_mask_array)
        """

        # We calculate the boolean mask arrays or check its validity if it was already given
        if mask_type=="raw":
            if(mask_ndarray.shape!=(self._Nz,self._Ny,self._Nx)):
                raise Exception("Wrong shape for the mask ndarray")
            mask = mask_ndarray
        else:
            if mask_type=="droplet":
                mask_formula = str((self._Lz/2)**2)+"-x**2-y**2-z**2"
            elif mask_type=="formula":
                mask_formula = mask_formula.replace("^","**")
            else:
                raise Exception(
                    "Unrecognised mask type: should be 'droplet', 'formula' or 'raw'")
            z,y,x = np.meshgrid(
                np.linspace(-0.5*self._dz*(self._Nz-1),0.5*self._dz*(self._Nz-1),self._Nz),
                np.linspace(-0.5*self._dy*(self._Ny-1),0.5*self._dy*(self._Ny-1),self._Ny),
                np.linspace(-0.5*self._dx*(self._Nx-1),0.5*self._dx*(self._Nx-1),self._Nx),
                indexing="ij")
            mask = eval(mask_formula)>=0

        # We check intersection with other subdomain
        other_subdomain_mask = self._vals>0
        if np.logical_and(other_subdomain_mask,mask).any() is True:
            raise Exception(
                "The new subdomain mask intersect with a previously defined subdomain")

        # We set the new subdomain id and add corresponding empty entries in the refractive
        # index array
        self._vals[mask] = self._refractive_indices.shape[0]
        self._refractive_indices = np.concatenate(
            (self._refractive_indices,np.array([["","",""]],dtype=object)),axis=0)


    def set_refractive_indices(self, \
            *, domain_id=0, niso=None, ne=None, no=None, n1=None, n2=None, n3=None):
        """Set the refractive indices of the default domain (domain_id=0) or one of the
        subdomain already specified (domain_id>0). In nemaktis, a refractive index is
        specified either as a simple constant float value (non-dispersive material) or a
        math string which can depend on the wavelength variable "lambda" (dispersive
        material, "lambda" is assumed to be in µm). 

        Depending on the phase symmetry of the subdomain, you should specify one or more
        refractive indices:

        .. code-block:: python

            # First case: isotropic medium -> one index of refraction
            sample.set_refractive_index(domain_id=..., niso=...)

            # Second case: uniaxial medium -> two indices of refraction
            sample.set_refractive_index(domain_id=..., ne=..., no=...)

            # Third case: biaxial medium -> three indices of refraction
            sample.set_refractive_index(domain_id=..., n1=..., n2=..., n3=...)

        For unixial media, the optical axis field "nfield" correspond to the extraordinary
        index ne. For biaxial media, the two optical axis fields "nfield" and "mfield"
        respectively correspond to the indices n1 and n2. 
        """
        if domain_id>=self._refractive_indices.shape[0]:
            raise Exception(
                "Cannot set indices, domain_id=%d was not yet created" % (domain_id,))
        if niso is not None: 
            self._refractive_indices[domain_id,:] = str(niso)
        elif ne is not None and no is not None:
            if self._optical_axis_field_type is None:
                raise Exception(
                    "You cannot specify unixial media if optical_axis_field_type is None")
            self._refractive_indices[domain_id,0] = str(ne)
            self._refractive_indices[domain_id,1:3] = str(no)
        elif n1 is not None and n2 is not None and n3 is not None:
            if self._optical_axis_field_type in [None,"Director","QTensor"]:
                raise Exception(
                    "Biaxiality is supported only if optical_axis_field_type is "
                    "DualDirector")
            self._refractive_indices[domain_id,0] = str(n1)
            self._refractive_indices[domain_id,1] = str(n2)
            self._refractive_indices[domain_id,2] = str(n3)


    def add_lower_isotropic_layer(self, *, nlayer, thickness):
        """
        Add an isotropic layer below the sample. Light is assumed to propagate in the
        z-direction, and will cross first the lower isotropic layers specified with this
        function, then the central regular cartesian mesh representing the sample, and
        finally the upper isotropic layers.

        The isotropic layers are always ordered in the direction of increasing z, so when
        you call this method, the new isotropic layer is inserted between the sample and all
        previously specified lower isotropic layers.

        Parameters
        ----------
        nlayer : float
            Refractive index of the new lower isotropic layer
        thickness : float
            Thickness (µm) of the new lower isotropic layer
        """
        self._lo_iso_layer_indices.append(nlayer)
        self._lo_iso_layer_thicknesses.append(thickness)


    def add_upper_isotropic_layer(self, *, nlayer, thickness):
        """
        Add an isotropic layer above the sample. Light is assumed to propagate in the
        z-direction, and will cross first the lower isotropic layers, then the central
        regular cartesian mesh representing the sample, and finally the upper isotropic
        layers specified with this function.

        The isotropic layers are always ordered in the direction of increasing z, so when
        you call this method, the new isotropic layer is inserted above all previously
        specified upper isotropic layers.

        Parameters
        ----------
        nlayer : float
            Refractive index of the new upper isotropic layer
        thickness : float
            Thickness (µm) of the new upper isotropic layer
        """
        self._up_iso_layer_indices.append(nlayer)
        self._up_iso_layer_thicknesses.append(thickness)


    def extend(self, scale_x, scale_y):
        """
        Extend the computational mesh in the ``xy`` plane by padding new points near the
        ``x`` and ``y`` boundaries. The associated new data points for all tensor fields
        included in this object (domain id field + all relevant optical axis fields) are
        initialized with the edge value of the tensor fields on the ``x`` and ``y``
        boundaries.

        Parameters
        ----------
        scale_x : float
            The mesh length in the x-direction will be scaled by this factor.
        scale_y : float
            The mesh length in the y-direction will be scaled by this factor.
        """
        super().extend(scale_x,scale_y)
        if self._optical_axis_field_type=="Director":
            self._nfield.extend(scale_x,scale_y)
        elif self._optical_axis_field_type=="DualDirector":
            self._nfield.extend(scale_x,scale_y)
            self._mfield.extend(scale_x,scale_y)
        elif self._optical_axis_field_type=="QTensor":
            self._qfield.extend(scale_x,scale_y)


    def rotate_90deg(self, axis):
        """
        Apply to all tensor fields included in this object (domain id field + all relevant
        optical axis fields) a solid rotation of 90 degrees around the specified axis. This
        is a lossless operation which does not rely on interpolation.
        
        Parameters
        ----------
        axis : str
            Axis around which to perform the rotation. Need to be under the form '[s]A'
            where the optional parameter 's'='+' or '-' decribes the sign of rotation and
            'A'='x', 'y' or 'z' defines the rotation axis.
        """
        super().rotate_90deg(axis)
        if self._optical_axis_field_type=="Director":
            self._nfield.rotate_90deg(axis)
        elif self._optical_axis_field_type=="DualDirector":
            self._nfield.rotate_90deg(axis)
            self._mfield.rotate_90deg(axis)
        elif self._optical_axis_field_type=="QTensor":
            self._qfield.rotate_90deg(axis)


    def rotate_180deg(self, axis):
        """
        Apply to all tensor fields included in this object (domain id field + all relevant
        optical axis fields) a solid rotation of 180 degrees around the specified axis. This
        is a lossless operation which does not rely on interpolation.
        
        Parameters
        ----------
        axis : str
            Axis around which to perform the rotation. Need to be under the form 'A' where
            'A'='x', 'y' or 'z' defines the rotation axis.
        """
        super().rotate_180deg(axis)
        if self._optical_axis_field_type=="Director":
            self._nfield.rotate_180deg(axis)
        elif self._optical_axis_field_type=="DualDirector":
            self._nfield.rotate_180deg(axis)
            self._mfield.rotate_180deg(axis)
        elif self._optical_axis_field_type=="QTensor":
            self._qfield.rotate_180deg(axis)


    def rotate(self, axis, angle, fill_value=None):
        """
        Apply to all tensor fields included in this object (domain id field + all relevant
        optical axis fields) a solid rotation of an arbitrary angle around the
        specified axis. This is a lossy operation that will rely on interpolation, so
        possible artefacts can appear if the tensor field data is not smooth enough.
        
        Note that the computational box stays unchanged after calling this method, which
        means that all mesh points falling out of the computational box after rotation will
        be clipped.  Please be aware that by default, new tensor values near the mesh
        boundaries will possibly be interpolated from extrapolated values outside the mesh
        since a regular cartesian mesh is not rotationally invariant. These values are
        generally spurious, and should always be ignored through the use of subdomain mask
        fully included inside the cartesian mesh whatever the applied rotation (e.g. a
        sphere centered inside the mesh). Alternatively, you can specify the fill_value
        parameter, which should contain a single array or a tuple of arrays (depending on
        the number optical axis fields contained in the object) for the uniform value(s) of
        the optical axis field(s) outside the cartesian mesh.
        
        Parameters
        ----------
        axis : str
            Axis around which to perform the rotation. Need to be under the form 'A' where
            'A'='x', 'y' or 'z' defines the rotation axis.
        angle : float
            Angle of rotation in degrees.
        fill_value : array or tuple of arrays, optional named parameter
            Value(s) of the optical axis field(s) outside the cartesian mesh.
        """
        super().rotate(axis,angle,method="nearest")
        if self._optical_axis_field_type=="Director":
            self._nfield.rotate(axis,angle,fill_value=fill_value)
        elif self._optical_axis_field_type=="DualDirector":
            if fill_value is None:
                fill_value = (None,None)
            self._nfield.rotate(axis,angle,fill_value=fill_value[0])
            self._mfield.rotate(axis,angle,fill_value=fill_value[1])
        elif self._optical_axis_field_type=="QTensor":
            self._qfield.rotate(axis,angle,fill_value=fill_value)


    def rescale_mesh(self, scaling_factor):
        """
        Uniformly scale the mesh of all tensor fields included in this object (domain id
        field + all relevant optical axis fields) using the given scaling factor.
        
        Parameters
        ----------
        scaling_factor : factor
            The mesh lengths and spacings will be multiplied by this factor.
        """
        super().rescale_mesh(scaling_factor)
        if self._optical_axis_field_type=="Director":
            self._nfield.rescale_mesh(scaling_factor)
        elif self._optical_axis_field_type=="DualDirector":
            self._nfield.rescale_mesh(scaling_factor)
            self._mfield.rescale_mesh(scaling_factor)
        elif self._optical_axis_field_type=="QTensor":
            self._qfield.rescale_mesh(scaling_factor)


    @property
    def nfield(self):
        """DirectorField object corresponding to the principal optical axis. This attribute
        exists only if the MicroscopeSample was initialized in Director or DualDirector
        mode."""
        if self._optical_axis_field_type in ["Director" or "DualDirector"]:
            return self._nfield
        else:
            raise Exception(
                "Cannot access nfield (optical_axis_field_type is %s)" \
                % (self._optical_axis_field_type))

    @nfield.setter
    def nfield(self, new_nfield):
        if self._optical_axis_field_type in ["Director" or "DualDirector"]:
            if not isinstance(new_nfield,DirectorField):
                raise Exception(
                    "The nfield attribute should be given a DirectorField object")
            elif self._vals.get_mesh_dimensions()!=new_nfield.get_mesh_dimensions \
                    or self.get_mesh_lengths()!=new_nfield.get_mesh_lengths():
                raise Exception(
                    "Incompatible mesh for the new nfield attribute")
            else:
                self._nfield = new_nfield
        else:
            raise Exception(
                "Cannot access nfield (optical_axis_field_type is %s)" \
                % (self._optical_axis_field_type))


    @property
    def mfield(self):
        """DirectorField object corresponding to the second optical axis. This attribute
        exists only if the MicroscopeSample was initialized in DualDirector mode."""
        if self._optical_axis_field_type=="DualDirector":
            return self._mfield
        else:
            raise Exception(
                "Cannot access mfield (optical_axis_field_type is %s)" \
                % (self._optical_axis_field_type))

    @mfield.setter
    def mfield(self, new_mfield):
        if self._optical_axis_field_type=="DualDirector":
            if not isinstance(new_mfield,DirectorField):
                raise Exception(
                    "The mfield attribute should be given a DirectorField object")
            elif self._vals.get_mesh_dimensions()!=new_mfield.get_mesh_dimensions \
                    or self.get_mesh_lengths()!=new_mfield.get_mesh_lengths():
                raise Exception(
                    "Incompatible mesh for the new mfield attribute")
            else:
                self._mfield = new_mfield
        else:
            raise Exception(
                "Cannot access mfield (optical_axis_field_type is %s)" \
                % (self._optical_axis_field_type))


    @property
    def qfield(self):
        """QtensorField object allowing to reconstruct the permittivity tensor. This
        attribute exists only if the MicroscopeSample was initialized in QTensor mode.""" 
        if self._optical_axis_field_type=="QTensor":
            return self._qfield
        else:
            raise Exception(
                "Cannot access qfield (optical_axis_field_type is %s)" \
                % (self._optical_axis_field_type))

    @qfield.setter
    def qfield(self, new_qfield):
        if self._optical_axis_field_type=="QTensor":
            if not isinstance(new_qfield,QTensorField):
                raise Exception(
                    "The qfield attribute should be given a QTensorField object")
            elif self._vals.get_mesh_dimensions()!=new_qfield.get_mesh_dimensions \
                    or self.get_mesh_lengths()!=new_qfield.get_mesh_lengths():
                raise Exception(
                    "Incompatible mesh for the new qfield attribute")
            else:
                self._qfield = new_qfield
        else:
            raise Exception(
                "Cannot access qfield (optical_axis_field_type is %s)" \
                % (self._optical_axis_field_type))


    def save_to_vti(self, file_name):
        """Save all the information stored in this object (domain id field, optical axis
        fields, refractive indices, isotropic layer thicknesses and indices) inside a vti
        file.

        Parameters
        ----------
        file_name : string 
            Path to the exported vti file. The ".vti" extension is automatically appended,
            no need to include it in this parameter (but in case you do only one extension
            will be added).
        """ 

        if file_name[-4:]==".vti":
            path = file_name
        else:
            path = file_name+".vti"
        print("{ Saving MicroscopeSample object to %s }" % path)

        vti_data = vtkImageData()
        vti_data.SetDimensions(self._Nx, self._Ny, self._Nz)
        vti_data.SetOrigin(-self._Lx/2, -self._Ly/2, -self._Lz/2)
        vti_data.SetSpacing(self._dx, self._dy, self._dz)

        tensor_data = \
            vn.numpy_to_vtk(self._vals.reshape((self._Nx*self._Ny*self._Nz,self._Nv)))
        tensor_data.SetName(self._name)
        vti_data.GetPointData().AddArray(tensor_data)

        if self._optical_axis_field_type=="Director":
            fields = [self._nfield]
        elif self._optical_axis_field_type=="DualDirector":
            fields = [self._nfield,self._mfield]
        elif self._optical_axis_field_type=="QTensor":
            fields = [self._qfield]

        for f in fields:
            tensor_data = \
                vn.numpy_to_vtk(f._vals.reshape((f._Nx*f._Ny*f._Nz,f._Nv)))
            tensor_data.SetName(f._name)
            vti_data.GetPointData().AddArray(tensor_data)

        domain_indices_data = vtkStringArray()
        for index_str in self._refractive_indices.flatten():
            if index_str=="":
                raise Exception(
                    "Invalid MicroscopeSample object: a refractive index was not set")
            domain_indices_data.InsertNextValue(index_str)
        domain_indices_data.SetName("domain_refractive_indices")
        vti_data.GetFieldData().AddArray(domain_indices_data)

        up_layer_indices_data = vn.numpy_to_vtk(self._up_iso_layer_indices)
        up_layer_indices_data.SetName("up_iso_layer_refractive_indices")
        vti_data.GetFieldData().AddArray(up_layer_indices_data)

        up_layer_thickness_data = vn.numpy_to_vtk(self._up_iso_layer_thicknesses)
        up_layer_thickness_data.SetName("up_iso_layer_thicknesses")
        vti_data.GetFieldData().AddArray(up_layer_thickness_data)

        lo_layer_indices_data = vn.numpy_to_vtk(self._lo_iso_layer_indices)
        lo_layer_indices_data.SetName("lo_iso_layer_refractive_indices")
        vti_data.GetFieldData().AddArray(lo_layer_indices_data)

        lo_layer_thickness_data = vn.numpy_to_vtk(self._lo_iso_layer_thicknesses)
        lo_layer_thickness_data.SetName("lo_iso_layer_thicknesses")
        vti_data.GetFieldData().AddArray(lo_layer_thickness_data)

        writer = vtkXMLImageDataWriter()
        writer.SetFileName(path)
        writer.SetInputData(vti_data)
        writer.Write()
