import sys
import numpy as np

from vtk import vtkXMLImageDataReader, vtkXMLImageDataWriter, vtkImageData
from vtk.util import numpy_support as vn
from scipy.interpolate import RectBivariateSpline, interp1d


# Function which allows to propagate exactly optical fields though a uniform birefringent
# slab of permittivity eps, provided their Fourier decomposion on the specified mesh is
# exact. Since the medium is homogeneous, we restrict the calculation to forward propagating
# modes.
def propagate_fields(input_field, mesh_lengths, mesh_dims, wavelength, eps):

    ################################
    # Mesh and fourier frequencies #
    ################################
    
    k0 = 2*np.pi/wavelength
    Lx,Ly,Lz = tuple(mesh_lengths)
    Nx,Ny,Nz = tuple(mesh_dims)
    
    xvals = np.linspace(-Lx/2,Lx/2,Nx)
    yvals = np.linspace(-Ly/2,Ly/2,Ny)
    zvals = np.linspace(0,Lz,Nz)
    
    dx = Lx/(Nx-1) if Nx!=1 else Lx
    dy = Ly/(Ny-1) if Ny!=1 else Ly
    dz = Lz/(Nz-1) if Nz!=1 else Lz
    
    px_vals = 2*np.pi*np.tile(np.fft.fftfreq(Nx, k0*dx)[np.newaxis,:],reps=(Ny,1))
    py_vals = 2*np.pi*np.tile(np.fft.fftfreq(Ny, k0*dy)[:,np.newaxis],reps=(1,Nx))
    p_vals = np.sqrt(px_vals**2+py_vals**2)
    
    
    #############################
    # Beam propagation matrices #
    #############################
    
    if np.allclose(eps,eps.T,atol=1e-8) is not True:
        print("Error: the permittivity tensor eps must be symmetric")
        sys.exit()
    
    # Transverse permittivity tensor
    eps_t = np.zeros((Ny,Nx,2,2))
    eps_t[:,:,0,0] = eps[0,0]-eps[0,2]**2/(eps[2,2]-p_vals**2)
    eps_t[:,:,0,1] = eps[0,1]-eps[0,2]*eps[1,2]/(eps[2,2]-p_vals**2)
    eps_t[:,:,1,0] = eps[0,1]-eps[0,2]*eps[1,2]/(eps[2,2]-p_vals**2)
    eps_t[:,:,1,1] = eps[1,1]-eps[1,2]**2/(eps[2,2]-p_vals**2)
    
    # NOTE: The exact wave eq in Fourier space can be defined as [d/dz^2 + I*R*d/dz + Q]E_perp=0,
    # with the operators Q and R defined below.
    
    # Phase+Diffraction operator
    Q_op = np.zeros((Ny,Nx,2,2))
    Q_op[:,:,0,0] = eps_t[:,:,0,0]-px_vals*(px_vals*eps_t[:,:,0,0]+py_vals*eps_t[:,:,1,0])/eps[2,2]-py_vals**2
    Q_op[:,:,0,1] = eps_t[:,:,0,1]-px_vals*(px_vals*eps_t[:,:,0,1]+py_vals*eps_t[:,:,1,1])/eps[2,2]+px_vals*py_vals
    Q_op[:,:,1,0] = eps_t[:,:,0,1]-py_vals*(px_vals*eps_t[:,:,0,0]+py_vals*eps_t[:,:,1,0])/eps[2,2]+px_vals*py_vals
    Q_op[:,:,1,1] = eps_t[:,:,1,1]-py_vals*(px_vals*eps_t[:,:,0,1]+py_vals*eps_t[:,:,1,1])/eps[2,2]-px_vals**2
    
    # Walk-off operator
    R_op = np.zeros((Ny,Nx,2,2))
    R_op[:,:,0,0] = 2*px_vals*eps[0,2]/eps[2,2] + \
                    px_vals*(eps[0,2]*py_vals**2-eps[1,2]*px_vals*py_vals)/(eps[2,2]*(eps[2,2]-p_vals**2))
    R_op[:,:,0,1] = (px_vals*eps[1,2]+py_vals*eps[0,2])/eps[2,2] + \
                    py_vals*(eps[0,2]*py_vals**2-eps[1,2]*px_vals*py_vals)/(eps[2,2]*(eps[2,2]-p_vals**2))
    R_op[:,:,1,0] = (px_vals*eps[1,2]+py_vals*eps[0,2])/eps[2,2] + \
                    px_vals*(eps[1,2]*px_vals**2-eps[0,2]*px_vals*py_vals)/(eps[2,2]*(eps[2,2]-p_vals**2))
    R_op[:,:,1,1] = 2*py_vals*eps[1,2]/eps[2,2] + \
                    py_vals*(eps[1,2]*px_vals**2-eps[0,2]*px_vals*py_vals)/(eps[2,2]*(eps[2,2]-p_vals**2))
    
    
    ##################################################
    # Eigenmode decomposition and evolution_operator #
    ##################################################
    
    # Light eigenmodes are computed from the operator [[0,nref I],[Q/nref, -R]] (this can be
    # shown by transforming the second-order wave equation with 2 dofs to a linear wave equation
    # with 4 dofs).
    
    # Reference index used to get a "good" scaling
    nref = np.sqrt(np.mean(np.linalg.eigvals(eps)))
    mat = np.zeros((Ny,Nx,4,4))
    mat[:,:,0:2,2:4] = nref*np.identity(2)
    mat[:,:,2:4,0:2] = Q_op/nref
    mat[:,:,2:4,2:4] = -R_op
    
    # We compute the eigenvalues and eigenvectors of mat and calculate the indices of forward modes
    eig_vals,eig_vecs = np.linalg.eig(mat)
    sorted_ids = np.argsort(np.real(eig_vals)+np.imag(eig_vals),axis=-1)[:,:,[-1,-2]]
    
    # We calculate the eigenvalues pz for the forward mode, and check that they fulfill the
    # stability criteria
    global_ids = (np.arange(0,Nx*Ny)[:,np.newaxis]*4+sorted_ids.reshape((Nx*Ny,2))).flatten()
    forward_pz = eig_vals.flatten()[global_ids].reshape(Ny,Nx,2)
    if np.any(np.imag(forward_pz)<-1e-8):
        print("Failure of selection algorithm: found an eigenvalue with negative imaginary part")
        sys.exit()
    
    # We calculate the basis-change matrices and their inverses from the components of the
    # eigenvectors
    p_ids = np.arange(0,Nx*Ny)[:,np.newaxis,np.newaxis]
    irow_ids = np.arange(0,2)[np.newaxis,:,np.newaxis]
    global_ids = (p_ids*16+irow_ids*4+sorted_ids.reshape((Nx*Ny,1,2))).flatten()
    trans_op = eig_vecs.flatten()[global_ids].reshape(Ny,Nx,2,2)
    inv_trans_op = np.linalg.inv(trans_op)
    
    # Finally, we calculate the evolution operator from the eigendecomposition operators
    phase_op = np.zeros((Ny,Nx,2,2),dtype=complex)
    phase_op[:,:,0,0] = np.exp(1j*forward_pz[:,:,0]*k0*dz)
    phase_op[:,:,1,1] = np.exp(1j*forward_pz[:,:,1]*k0*dz)
    evolution_op = trans_op @ phase_op @ inv_trans_op
    
    
    #####################
    # Field propagation #
    #####################

    fft_field = np.zeros((Nz,Ny,Nx,2),dtype=complex)
    fft_field[0,:,:,:] = np.fft.fft2(input_field,axes=(0,1))
    for iz in range(1,Nz):
        fft_field[iz,:,:,:] = np.squeeze(
            evolution_op @ fft_field[iz-1,:,:,:,np.newaxis], axis=-1)
    field = np.fft.ifft2(fft_field,axes=(1,2))

    return field


# Allows to import an input profile from a bulk simulation file exported by Nemaktis.
# The parameter pol="X" or "Y" specifies the input polarisation of the beam. 
def import_input_field(vti_file, mesh_lengths, mesh_dims,*, pol):
    if pol is not "X" and pol is not "Y":
        print("Error: the parameter pol should be \"X\" or \"Y\"")
        return

    reader = vtkXMLImageDataReader()
    reader.SetFileName(vti_file)
    reader.Update()
    
    point_data = reader.GetOutput().GetPointData()
    dims = np.array(reader.GetOutput().GetDimensions())
    origin = np.array(reader.GetOutput().GetOrigin())
    spacings = np.array(reader.GetOutput().GetSpacing())

    field_dim = point_data.GetArray("E_real_inputX_0_0").GetNumberOfComponents()
    
    field_vals = \
        vn.vtk_to_numpy(point_data.GetArray("E_real_input%s_0_0" % (pol,))) + \
        1j*vn.vtk_to_numpy(point_data.GetArray("E_imag_input%s_0_0" % (pol,)))
    Ex_input_vals = field_vals.reshape((dims[2],dims[1],dims[0],field_dim))[0,:,:,0]
    Ey_input_vals = field_vals.reshape((dims[2],dims[1],dims[0],field_dim))[0,:,:,1]

    if dims[0]>3 and dims[1]>3:
        xvals = np.linspace(origin[0],origin[0]+(dims[0]-1)*spacings[0],dims[0])
        yvals = np.linspace(origin[1],origin[1]+(dims[1]-1)*spacings[1],dims[1])

        Ex_input_real = RectBivariateSpline(yvals,xvals,np.real(Ex_input_vals))
        Ex_input_imag = RectBivariateSpline(yvals,xvals,np.imag(Ex_input_vals))
        Ey_input_real = RectBivariateSpline(yvals,xvals,np.real(Ey_input_vals))
        Ey_input_imag = RectBivariateSpline(yvals,xvals,np.imag(Ey_input_vals))

        new_xvals = np.linspace(-mesh_lengths[0]/2,mesh_lengths[0]/2,mesh_dims[0])
        new_yvals = np.linspace(-mesh_lengths[1]/2,mesh_lengths[1]/2,mesh_dims[1])

        input_field = np.zeros((mesh_dims[1],mesh_dims[0],2),dtype=complex)
        input_field[:,:,0] = Ex_input_real(new_yvals,new_xvals)+1j*Ex_input_imag(new_yvals,new_xvals)
        input_field[:,:,1] = Ey_input_real(new_yvals,new_xvals)+1j*Ey_input_imag(new_yvals,new_xvals)
    
    elif dims[0]>3:
        xvals = np.linspace(origin[0],origin[0]+(dims[0]-1)*spacings[0],dims[0])

        Ex_input_real = interp1d(xvals,np.real(Ex_input_vals[int(dims[1]/2),:]),kind="cubic")
        Ex_input_imag = interp1d(xvals,np.imag(Ex_input_vals[int(dims[1]/2),:]),kind="cubic")
        Ey_input_real = interp1d(xvals,np.real(Ey_input_vals[int(dims[1]/2),:]),kind="cubic")
        Ey_input_imag = interp1d(xvals,np.imag(Ey_input_vals[int(dims[1]/2),:]),kind="cubic")

        new_xvals = np.linspace(-mesh_lengths[0]/2,mesh_lengths[0]/2,mesh_dims[0])

        input_field = np.zeros((mesh_dims[1],mesh_dims[0],2),dtype=complex)
        input_field[:,:,0] = (Ex_input_real(new_xvals)+1j*Ex_input_imag(new_xvals))[np.newaxis,:]
        input_field[:,:,1] = (Ey_input_real(new_xvals)+1j*Ey_input_imag(new_xvals))[np.newaxis,:]
    
    elif dims[1]>3:
        yvals = np.linspace(origin[1],origin[1]+(dims[1]-1)*spacings[1],dims[1])

        Ex_input_real = interp1d(yvals,np.real(Ex_input_vals[:,int(dims[0]/2)]),kind="cubic")
        Ex_input_imag = interp1d(yvals,np.imag(Ex_input_vals[:,int(dims[0]/2)]),kind="cubic")
        Ey_input_real = interp1d(yvals,np.real(Ey_input_vals[:,int(dims[0]/2)]),kind="cubic")
        Ey_input_imag = interp1d(yvals,np.imag(Ey_input_vals[:,int(dims[0]/2)]),kind="cubic")

        new_yvals = np.linspace(-mesh_lengths[1]/2,mesh_lengths[1]/2,mesh_dims[1])

        input_field = np.zeros((mesh_dims[1],mesh_dims[0],2),dtype=complex)
        input_field[:,:,0] = (Ex_input_real(new_yvals)+1j*Ex_input_imag(new_yvals))[:,np.newaxis]
        input_field[:,:,1] = (Ey_input_real(new_yvals)+1j*Ey_input_imag(new_yvals))[:,np.newaxis]

    return input_field


# Allows to import the complete 3D optical fields from a bulk simulation file exported by
# Nemaktis. The parameter pol="X" or "Y" specifies the input polarisation of the beam. 
def import_field(vti_file, *, pol):
    if pol is not "X" and pol is not "Y":
        print("Error: the parameter pol should be \"X\" or \"Y\"")
        return

    reader = vtkXMLImageDataReader()
    reader.SetFileName(vti_file)
    reader.Update()
    
    point_data = reader.GetOutput().GetPointData()
    dims = np.array(reader.GetOutput().GetDimensions())
    origin = np.array(reader.GetOutput().GetOrigin())
    spacings = np.array(reader.GetOutput().GetSpacing())

    field_dim = point_data.GetArray("E_real_inputX_0_0").GetNumberOfComponents()
    field_vals = \
        vn.vtk_to_numpy(point_data.GetArray("E_real_input%s_0_0" % (pol,))) + \
        1j*vn.vtk_to_numpy(point_data.GetArray("E_imag_input%s_0_0" % (pol,)))
    return field_vals.reshape((dims[2],dims[1],dims[0],field_dim))[:,:,:,[0,1]]
