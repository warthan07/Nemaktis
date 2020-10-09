import time
import os
import json

import dtmm
dtmm.conf.set_fftlib("mkl_fft")

import numpy as np
import nemaktis as nm
import matplotlib.pyplot as plt

from copy import deepcopy
from propagate_fields import *
from scipy.interpolate import RegularGridInterpolator


########################
# Simulation constants #
########################
	
ne = 1.75
no = 1.5
wavelength = 0.6
k0 = 2*np.pi/wavelength

L = 10

theta = np.pi/4
phi = 0
def nx(x,y,z):
    return np.cos(phi)*np.sin(theta)*np.ones(x.shape)
def ny(x,y,z):
    return np.sin(phi)*np.sin(theta)*np.ones(x.shape)
def nz(x,y,z):
    return np.cos(theta)*np.ones(x.shape)

for diffraction in [1,3,5,7]:
    ####################
    # DTMM simulations #
    ####################
    
    N_vals = np.logspace(np.log10(11), np.log10(301), 20).astype(int)
    d_vals = L/(N_vals-1)
    elapsed_times = np.zeros(len(N_vals))
    
    if not os.path.exists('results'):
        os.makedirs('results')
    
    f = open("results/elapsed_times_dtmm_diff=%d.dat" % (diffraction,), "w")
    f.write("# N\td\tt\n")
    for iN in range(0,len(N_vals)):
        N = N_vals[iN]
    
        print("")
        print("##################")
        print("### N=%d ###" % (N,))
        print("##################")
    
        nfield = nm.DirectorField(
            mesh_lengths=(L, L, L), mesh_dimensions=(N, N, N))
        nfield.init_from_funcs(nx,ny,nz)
    
        dims = nfield.get_mesh_dimensions()
        spacings = nfield.get_mesh_spacings()
    
        optical_data = dtmm.director2data(
           nfield.vals, no = no, ne = ne, nhost = no,
           thickness = spacings[2]/spacings[1]*np.ones(dims[2]))
    
        t0 = time.perf_counter()
    
        gaussian_window = dtmm.window.gaussian_beam(
            (dims[1],dims[0]), 1/(np.sqrt(2)*spacings[1]), k0)
        field_data_in = dtmm.illumination_data(
            (dims[1],dims[0]), [wavelength], pixelsize=spacings[1], n=no,
            window = gaussian_window)
        field_data_out = dtmm.transfer_field(
            field_data_in, optical_data, nin=no, diffraction=diffraction, ret_bulk=True)[0]
    
        elapsed_times[iN] = time.perf_counter()-t0
        f.write("%d\t%f\t%f\n" % (N,d_vals[iN],elapsed_times[iN]))
    
        # In order to be consistent with the BPM backend, we absorb the phase factor k0*no*z
        z_vals = spacings[2]*np.arange(0,N+2)[:,np.newaxis,np.newaxis,np.newaxis,np.newaxis,np.newaxis]
        field_data_out *= np.exp(-1j*k0*no*z_vals)
    
        bulk_filename = "results/dtmm_sol_diff=%d_N=%d" % (diffraction,N,)
        print("{ Saving optical fields to "+bulk_filename+".vti }")
        lengths = nfield.get_mesh_lengths()
    
        vti_data = vtkImageData()
        vti_data.SetDimensions(dims[0], dims[1], dims[2])
        vti_data.SetOrigin(-lengths[0]/2, -lengths[1]/2, -lengths[2]/2)
        vti_data.SetSpacing(spacings[0], spacings[1], spacings[2])
    
        Np = dims[0]*dims[1]*dims[2]
        E_inputX = field_data_out[:-2,0,0,[0,2],:,:].transpose(
            (1,0,2,3)).reshape((2,Np)).transpose()
        E_inputY = field_data_out[:-2,1,0,[0,2],:,:].transpose(
            (1,0,2,3)).reshape((2,Np)).transpose()
    
        E_real_inputX = vn.numpy_to_vtk(np.real(E_inputX))
        E_real_inputX.SetName("E_real_inputX_0_0")
        vti_data.GetPointData().AddArray(E_real_inputX)
        E_imag_inputX = vn.numpy_to_vtk(np.imag(E_inputX))
        E_imag_inputX.SetName("E_imag_inputX_0_0")
        vti_data.GetPointData().AddArray(E_imag_inputX)
    
        E_real_inputY = vn.numpy_to_vtk(np.real(E_inputY))
        E_real_inputY.SetName("E_real_inputY_0_0")
        vti_data.GetPointData().AddArray(E_real_inputY)
        E_imag_inputY = vn.numpy_to_vtk(np.imag(E_inputY))
        E_imag_inputY.SetName("E_imag_inputY_0_0")
        vti_data.GetPointData().AddArray(E_imag_inputY)
    
        writer = vtkXMLImageDataWriter()
        writer.SetFileName(bulk_filename+".vti")
        writer.SetInputData(vti_data)
        writer.Write()
    
    f.close()
    
        
    ################################################
    # Calculation of exact solution based on FFT #
    ################################################
    
    print("")
    print("#################################")
    print("# Calculation of exact solution #")
    print("#################################")
    
    optical_axis = np.array([np.cos(phi)*np.sin(theta),np.sin(phi)*np.sin(theta),np.cos(theta)]).reshape((3,1))
    eps = no**2*np.identity(3)+(ne**2-no**2)*np.kron(optical_axis,optical_axis.T)
    
    N = N_vals[-1]
    mesh_lengths = [L,L,L]
    mesh_dims = [N,N,N]
    
    input_field = import_input_field(
        "results/dtmm_sol_diff=%d_N=%d.vti" % (diffraction,N_vals[-1],), mesh_lengths, mesh_dims, pol="X")
    exact_field_vals = propagate_fields(
        input_field, mesh_lengths, mesh_dims, wavelength, eps)
    
    z_vals = np.linspace(0,L,N)[:,np.newaxis,np.newaxis,np.newaxis]
    exact_field_vals = np.squeeze(exact_field_vals*np.exp(-1j*k0*no*z_vals))
    fft_exact_field_vals = np.fft.fft2(exact_field_vals,axes=(1,2))
    fft_exact_ampl_vals = np.sqrt(np.sum(np.abs(fft_exact_field_vals)**2,axis=-1))
    
    xs_ref = np.linspace(-L/2,L,N)
    ys_ref = np.linspace(-L/2,L,N)
    zs_ref = np.linspace(0,L,N)
    Zs_ref,Ys_ref,Xs_ref = np.meshgrid(zs_ref, ys_ref, xs_ref, indexing="ij")
    pts_ref = np.stack((Zs_ref,Ys_ref,Xs_ref),axis=-1)
    
    qx_vals = np.tile(2*np.pi*np.fft.fftfreq(N,L/(N-1))/k0, (N,N,1))
    qy_vals = np.tile(2*np.pi*np.fft.fftfreq(N,L/(N-1))/k0, (N,N,1)).transpose((0,2,1))
    q_vals = np.sqrt(qx_vals**2+qy_vals**2)
    
    NA = 0.4
    mask = q_vals**2 < NA**2*(1-(q_vals/no)**2)
    
    
    ##########################################
    # Calculation of the computational error #
    ##########################################
    
    class complex_interp:
        def __init__(self, zs, ys, xs, vals):
            self.interp_real = RegularGridInterpolator((zs, ys, xs), np.real(vals))
            self.interp_imag = RegularGridInterpolator((zs, ys, xs), np.imag(vals))
        def __call__(self, points):
            return self.interp_real(points)+1j*self.interp_imag(points)
    
    print("")
    print("#######################################")
    print("# Calculation of computational errors #")
    print("#######################################")
    
    max_errors = np.zeros(elapsed_times.shape)
    min_errors = np.zeros(elapsed_times.shape)
    mean_errors = np.zeros(elapsed_times.shape)
    with open("results/errors_and_times_dtmm_3D_diff=%d.dat" % (diffraction,), mode='w') as f:
        f.write("# N\td\tt\terr_mean\terr_min\terr_max\n")
        for iN in range(0,len(N_vals)):
            N = N_vals[iN]
    
            field_vals = import_field("results/dtmm_sol_diff=%d_N=%d.vti" % (diffraction,N,), pol="X")
            xs = np.linspace(-L/2,L,N)
            ys = np.linspace(-L/2,L,N)
            zs = np.linspace(0,L,N)
    
            Ex = complex_interp(zs, ys, xs, field_vals[:,:,:,0])
            Ey = complex_interp(zs, ys, xs, field_vals[:,:,:,1])
            fft_Ex_err = np.fft.fft2(exact_field_vals[:,:,:,0]-Ex(pts_ref),axes=(1,2))
            fft_Ey_err = np.fft.fft2(exact_field_vals[:,:,:,1]-Ey(pts_ref),axes=(1,2))
            fft_err = \
                np.sqrt(np.abs(fft_Ex_err)**2+np.abs(fft_Ey_err)**2) / \
                (1e-8+2*np.abs(fft_exact_ampl_vals))

            max_errors[iN] = np.max(fft_err[mask])
            min_errors[iN] = np.min(fft_err[mask])
            mean_errors[iN] = np.mean(fft_err[mask])
            f.write("%d\t%f\t%f\t%f\t%f\t%f\n" % (
                N,d_vals[iN],elapsed_times[iN],mean_errors[iN],min_errors[iN],max_errors[iN]))
    
    #  plt.loglog(d_vals,mean_errors)
    #  plt.xlabel("Mesh spacing dy (µm)")
    #  plt.ylabel("Computational error for DTMM backend")

    #  plt.show()
