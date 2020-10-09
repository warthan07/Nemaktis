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


########################
# Simulation constants #
########################
	
ne = 1.75
no = 1.5
wavelength = 0.6
k0 = 2*np.pi/wavelength

Ly = 10
Lz = 10

angle = np.pi/4
def nx(x,y,z):
    return np.zeros(x.shape)
def ny(x,y,z):
    return -np.sin(angle)*np.ones(x.shape)
def nz(x,y,z):
    return np.cos(angle)*np.ones(x.shape)

for diffraction in [1,3,5,7]:
    ####################
    # DTMM simulations #
    ####################
    
    Ny_vals = np.logspace(np.log10(31), np.log10(1001), 20).astype(int)
    dy_vals = Ly/(Ny_vals-1)
    elapsed_times = np.zeros(len(Ny_vals))
    
    if not os.path.exists('results'):
        os.makedirs('results')
    
    f = open("results/elapsed_times_dtmm_diff=%d.dat" % (diffraction,), "w")
    f.write("# Ny\tdy\tt\n")
    for iN in range(0,len(Ny_vals)):
        Ny = Ny_vals[iN]
        Nz = Ny
    
        print("")
        print("##################")
        print("### Ny=%d ###" % (Ny,))
        print("##################")
    
        nfield = nm.DirectorField(
            mesh_lengths=(1, Ly, Lz), mesh_dimensions=(1, Ny, Nz))
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
        f.write("%d\t%f\t%f\n" % (Ny,dy_vals[iN],elapsed_times[iN]))
    
        # In order to be consistent with the BPM backend, we absorb the phase factor k0*no*z
        z_vals = spacings[2]*np.arange(0,Nz+2)[:,np.newaxis,np.newaxis,np.newaxis,np.newaxis,np.newaxis]
        field_data_out *= np.exp(-1j*k0*no*z_vals)
    
        bulk_filename = "results/dtmm_sol_diff=%d_Ny=%d" % (diffraction,Ny,)
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
    
    optical_axis = np.array([0,-np.sin(angle),np.cos(angle)]).reshape((3,1))
    eps = no**2*np.identity(3)+(ne**2-no**2)*np.kron(optical_axis,optical_axis.T)
    
    Ny = Ny_vals[-1]
    Nz = Ny
    mesh_lengths = [1,Ly,Lz]
    mesh_dims = [1,Ny,Nz]
    
    input_field = import_input_field(
        "results/dtmm_sol_diff=%d_Ny=%d.vti" % (diffraction,Ny_vals[-1],), mesh_lengths, mesh_dims, pol="Y")
    exact_field_vals = propagate_fields(
        input_field, mesh_lengths, mesh_dims, wavelength, eps)
    
    z_vals = np.linspace(0,Lz,Nz)[:,np.newaxis,np.newaxis,np.newaxis]
    exact_field_vals = np.squeeze(exact_field_vals*np.exp(-1j*k0*no*z_vals))
    fft_exact_field_vals = np.fft.fft(exact_field_vals,axis=1)
    fft_exact_ampl_vals = np.sqrt(np.sum(np.abs(fft_exact_field_vals)**2,axis=-1))
    
    ys_ref = np.linspace(-Ly/2,Ly,Ny)
    zs_ref = np.linspace(0,Lz,Nz)
    q_vals = np.tile(2*np.pi*np.fft.fftfreq(Ny,Ly/(Ny-1))/k0, (Nz,1))
    
    NA = 0.4
    mask = q_vals**2 < NA**2*(1-(q_vals/no)**2)
    
    
    ##########################################
    # Calculation of the computational error #
    ##########################################
    
    class complex_interp:
        def __init__(self, zs, ys, vals, kx, ky):
            self.interp_real = RectBivariateSpline(zs, ys, np.real(vals), kx=kx, ky=ky)
            self.interp_imag = RectBivariateSpline(zs, ys, np.imag(vals), kx=kx, ky=ky)
        def __call__(self, zs, ys):
            return self.interp_real(zs,ys)+1j*self.interp_imag(zs,ys)
    
    print("")
    print("#######################################")
    print("# Calculation of computational errors #")
    print("#######################################")
    
    max_errors = np.zeros(elapsed_times.shape)
    min_errors = np.zeros(elapsed_times.shape)
    mean_errors = np.zeros(elapsed_times.shape)
    with open("results/errors_and_times_dtmm_2D_diff=%d.dat" % (diffraction,), mode='w') as f:
        f.write("# Ny\tdy\tt\terr_mean\terr_min\terr_max\n")
        for iN in range(0,len(Ny_vals)):
            Ny = Ny_vals[iN]
            Nz = Ny
    
            field_vals = import_field("results/dtmm_sol_diff=%d_Ny=%d.vti" % (diffraction,Ny,), pol="Y")[:,:,0,:]
            ys = np.linspace(-Ly/2,Ly,Ny)
            zs = np.linspace(0,Lz,Nz)
    
            Ex = complex_interp(zs, ys, field_vals[:,:,0], kx=3, ky=3)
            Ey = complex_interp(zs, ys, field_vals[:,:,1], kx=3, ky=3)
            fft_Ex_err = np.fft.fft(exact_field_vals[:,:,0]-Ex(zs_ref,ys_ref),axis=-1)
            fft_Ey_err = np.fft.fft(exact_field_vals[:,:,1]-Ey(zs_ref,ys_ref),axis=-1)
            fft_err = \
                np.sqrt(np.abs(fft_Ex_err)**2+np.abs(fft_Ey_err)**2) / \
                (1e-8+2*np.abs(fft_exact_ampl_vals))
    
            max_errors[iN] = np.max(fft_err[mask])
            min_errors[iN] = np.min(fft_err[mask])
            mean_errors[iN] = np.mean(fft_err[mask])
            f.write("%d\t%f\t%f\t%f\t%f\t%f\n" % (
                Ny,dy_vals[iN],elapsed_times[iN],mean_errors[iN],min_errors[iN],max_errors[iN]))
    
    #  plt.subplot(1,2,1)
    #  plt.loglog(dy_vals,mean_errors)
    #  plt.xlabel("Mesh spacing dy (µm)")
    #  plt.ylabel("Computational error for DTMM backend")

    #  plt.subplot(1,2,2)
    #  plt.imshow(np.abs(exact_field_vals[:,:,1]),extent=(-Ly/2,Ly/2,0,Lz),origin="lower")
    #  plt.title("Amplitude of Ey")
    #  plt.colorbar()

    #  plt.show()
