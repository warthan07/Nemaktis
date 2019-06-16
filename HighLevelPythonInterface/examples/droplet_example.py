import nemaktis as nm
import numpy as np
import os.path

if os.path.isfile("optical_fields.vti"):
    # If the optical field were already calculated and exported by this
    # script, we directly load them
    output_fields = nm.OpticalFields(vti_file="optical_fields.vti")

else:
    # Else, we need to load a director field and propagate fields
    # through it. We use a simple ansatz for a double twist droplet.
    q = 2*np.pi/20
    def nx(x,y,z):
        r = np.sqrt(x**2+y**2)
        return -q*y*np.sinc(q*r)
    def ny(x,y,z):
        r = np.sqrt(x**2+y**2)
        return q*x*np.sinc(q*r)
    def nz(x,y,z):
        r = np.sqrt(x**2+y**2)
        return np.cos(q*r)

    nfield = nm.DirectorField(
        mesh_lengths=(10, 10, 10), mesh_dimensions=(80, 80, 80))
    nfield.init_from_funcs(nx,ny,nz)
    nfield.normalize()

    nfield.rotate_90deg("x")
    nfield.extend(2,2)
    nfield.set_mask(mask_type="droplet")

    # The LCMaterial object contains the details of the materials
    # of the LC sample: LC layer + possible isotropic layers above it
    # (a glass plate for example)
    mat = nm.LCMaterial(
        director_field=nfield,
        ne="1.6933+0.0078/lambda^2+0.0028/lambda^4",
        no="1.4990+0.0072/lambda^2+0.0003/lambda^4",
        nhost=1.55)
    mat.add_isotropic_layer(nlayer=1.51, thickness=1000)

    wavelengths = np.linspace(0.4, 0.8, 11)
    sim = nm.LightPropagator(material=mat, wavelengths=wavelengths, numerical_aperture=0.4)
    output_fields = sim.propagate_fields(method="dtmm")

    # We save the optical fields in a vti file
    output_fields.save_to_vti("optical_fields")

# Finally, the optical fields are visualized as in a real microscope
viewer = nm.FieldViewer(output_fields)
viewer.plot()
