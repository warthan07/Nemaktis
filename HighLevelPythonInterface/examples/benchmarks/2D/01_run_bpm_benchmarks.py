import time
import os
import json

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


###################
# BPM simulations #
###################

json_params_base = {
        "Algorithm settings": {
            "General": {
                "LC field type":        "Director",
                "Results folder name":  "results" },
            "Beam propagation": {
                "Number of substeps per slab":  1,
                "N Woodbury steps":             4 }},
        "Physics settings": {
            "Initial conditions": {
                "Beam profile":      "GaussianBeam(0.5)",
                "LC field file":     "nfield.vti",
                "Basis convention":  "XYZ" },
            "Coefficients": {
                "no":                str(no),
                "ne":                str(ne),
                "nhost":             str(no),
                "nin":               str(no),
                "Mean wavelength":   wavelength,
                "Spectral FWHM":     0.2,
                "N wavelengths":     1,
                "Condenser numerical aperture": 0,
                "N radial illumination directions": 1 }},
        "Postprocessor settings": {
            "Bulk output": {
                "Activate":   True,
                "Base name":  "" },
            "Screen output": {
                "Activate":                           False,
                "Base name":                          "screen",
                "Isotropic layer thicknesses":        [],
                "Isotropic layer refractive indices": [],
                "Focalisation z-shift":               0,
                "Numerical aperture":                 0.3 }}}

Ny_vals = np.logspace(np.log10(31), np.log10(1001), 20).astype(int)
dy_vals = Ly/(Ny_vals-1)
elapsed_times = np.zeros(len(Ny_vals))

if not os.path.exists('results'):
    os.makedirs('results')

f = open("results/elapsed_times_bpm.dat", "w")
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
    nfield.save_to_vti("nfield.vti")

    json_params = deepcopy(json_params_base)
    json_params["Postprocessor settings"]["Bulk output"]["Base name"] = \
        "bpm_sol_Ny=%d" % (Ny,)

    json_file = open("settings.par","w")
    json.dump(json_params, json_file)
    json_file.close()

    t0 = time.perf_counter()
    os.system("bpm-solver -x settings.par")
    elapsed_times[iN] = time.perf_counter()-t0

    f.write("%d\t%f\t%f\n" % (Ny,dy_vals[iN],elapsed_times[iN]))
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
    "results/bpm_sol_Ny=%d.vti" % (Ny,), mesh_lengths, mesh_dims, pol="Y")
exact_field_vals = propagate_fields(
    input_field, mesh_lengths, mesh_dims, wavelength, eps)

# With the BPM backend of Nemaktis, the optical fields are defined up to a xy-invariant
# phase term in k0*no*z, so we take into account this phase term here to allow an easy
# comparison between the two calculations.
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
with open("results/errors_and_times_bpm_2D.dat", mode='w') as f:
    f.write("# Ny\tdy\tt\terr_mean\terr_min\terr_max\n")
    for iN in range(0,len(Ny_vals)):
        Ny = Ny_vals[iN]
        Nz = Ny

        field_vals = np.squeeze(import_field("results/bpm_sol_Ny=%d.vti" % (Ny,), pol="Y"))
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

plt.subplot(1,2,1)
plt.loglog(dy_vals,mean_errors)
plt.xlabel("Mesh spacing dy (µm)")
plt.ylabel("Computational error for BPM backend")

plt.subplot(1,2,2)
plt.imshow(np.abs(exact_field_vals[:,:,1]),extent=(-Ly/2,Ly/2,0,Lz),origin="lower")
plt.title("Amplitude of Ey")
plt.colorbar()

plt.show()
