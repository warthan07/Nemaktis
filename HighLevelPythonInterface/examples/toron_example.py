import nemaktis as nm

import numpy as np
from numpy import pi, exp, cos, sin, sign, sqrt
import os.path

if os.path.isfile("director_field.vti"):
    # If the director fiel was already calculated and exported by this
    # script, we directly load it
    nfield = nm.DirectorField(vti_file="director_field.vti")

else:
    # Else, we calculate the director field of a toron from an ansatz
    R = 5
    d = 20
    xi = 0.5
    def alpha(r,z):
        u = 1 if abs(3*z/d)<1 else 0
        return pi/2*( (1-exp(-r**2/(2*xi**2)))*cos(pi*z/d) + u*exp(-r**2/(2*xi**2)) )
    def beta(z):
        return -pi*abs(z)/d
    def gamma(r):
        return pi*(1-exp(-r**2/R**2))
    def nr(r,z):
        a = alpha(r,z)
        b = beta(z)
        g = gamma(r)
        return sign(z)*(sin(a)*sin(b)*sin(g)-sin(2*a)*cos(b)*cos(g/2)**2)
    def nphi(r,z):
        a = alpha(r,z)
        b = beta(z)
        g = gamma(r)
        return sin(a)*cos(b)*sin(g)+sin(2*a)*sin(b)*cos(g/2)**2
    def nx(x,y,z):
        r = sqrt(x**2+y**2)
        return (x/r)*nr(r,z)-(y/r)*nphi(r,z) if r!=0 else 0
    def ny(x,y,z):
        r = sqrt(x**2+y**2)
        return (y/r)*nr(r,z)+(x/r)*nphi(r,z) if r!=0 else 0
    def nz(x,y,z):
        r = sqrt(x**2+y**2)
        a = alpha(r,z)
        g = gamma(r)
        return 1-2*sin(a)**2*cos(g/2)**2
    
    nfield = nm.DirectorField(
        mesh_lengths=(40, 40, 20), mesh_dimensions=(100, 100, 100))
    nfield.init_from_funcs(nx,ny,nz)
    nfield.save_to_vti("director_field")

    
# We propagate fields through the lc layer
mat = nm.LCMaterial(
    lc_field=nfield,
    ne=1.50,
    no=1.65,
    nhost=1.55)

wavelengths = np.linspace(0.4, 0.8, 11)
sim = nm.LightPropagator(
    material=mat, wavelengths=wavelengths, max_NA_objective=0.4)
output_fields = sim.propagate_fields(method="dtmm")

# Finally, the optical fields are visualized as in a real microscope
viewer = nm.FieldViewer(output_fields)
viewer.plot()
