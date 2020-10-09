import time
import os
import json

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

N_vals = np.logspace(np.log10(11), np.log10(301), 20).astype(int)
d_vals = L/(N_vals-1)
elapsed_times = np.zeros(len(N_vals))

if not os.path.exists('results'):
    os.makedirs('results')

f = open("results/elapsed_times_bpm.dat", "w")
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
    nfield.save_to_vti("nfield.vti")

    json_params = deepcopy(json_params_base)
    json_params["Postprocessor settings"]["Bulk output"]["Base name"] = \
        "bpm_sol_N=%d" % (N,)

    json_file = open("settings.par","w")
    json.dump(json_params, json_file)
    json_file.close()

    t0 = time.perf_counter()
    os.system("bpm-solver -x settings.par")
    elapsed_times[iN] = time.perf_counter()-t0

    f.write("%d\t%f\t%f\n" % (N,d_vals[iN],elapsed_times[iN]))
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
    "results/bpm_sol_N=%d.vti" % (N,), mesh_lengths, mesh_dims, pol="X")
exact_field_vals = propagate_fields(
    input_field, mesh_lengths, mesh_dims, wavelength, eps)

# With the BPM backend of Nemaktis, the optical fields are defined up to a xy-invariant
# phase term in k0*no*z, so we take into account this phase term here to allow an easy
# comparison between the two calculations.
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
with open("results/errors_and_times_bpm_3D.dat", mode='w') as f:
    f.write("# N\td\tt\terr_mean\terr_min\terr_max\n")
    for iN in range(0,len(N_vals)):
        N = N_vals[iN]

        field_vals = np.squeeze(import_field("results/bpm_sol_N=%d.vti" % (N,), pol="X"))
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

plt.loglog(d_vals,mean_errors)
plt.xlabel("Mesh spacing dy (µm)")
plt.ylabel("Computational error for BPM backend")

plt.show()
