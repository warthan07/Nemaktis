.. _tutorial:

Tutorial
========

This tutorial provides a hands-on introduction to the python package ``nemaktis``.
You will learn the different ways of creating director field data,
how to define the sample geometry and material constants, and how to propagate and
visualise optical fields.

First of all, open your favorite text/code editor and create a new python file
(which we will call ``script.py`` in the following). The script can be tested at any
moment in a terminal on condition that the conda environment in which you installed
``nemaktis`` is activated (``conda activate [environment name]``): ::

    cd [path to your script]
    python script.py

Alternatively, you can work interactively with ipython (which must be run from a terminal in
which the conda environment for ``nemaktis`` is activated).


.. _nfield:

Defining a DirectorField
-------------------------

Before starting using ``nemaktis``, we of course need to import the associated python package.
We will also import numpy, which will be needed to define arrays:

.. code-block:: python

    import nemaktis as nm
    import numpy as np

Next, we need to define a director field. In ``nemaktis``, any vector field is represented
internally on a cartesian regular mesh as a numpy array of shape ``(Nz,Ny,Nx,Nv)``, where
``Nv`` is the dimension of the vector data (3 for a director field, 6 for a symmétric tensor)
and ``Nx``, ``Ny`` and ``Nz`` are the number of mesh points in each spatial direction. In
addition to these variables, one needs to specify the total lengths of the mesh in each spatial
direction, which we will call ``Lx``, ``Ly`` and ``Lz`` in the following. All lengths are
in micrometer in ``nemaktis``, and the mesh for the director field is always centerered on the
origin (which means that the spatial coordinate ``u=x,y,z`` is always running from ``-Lu/2``
to ``Lu/2``).

Here, we will start by defining an empty director field object on a mesh of dimensions ``80x80x80``
and lengths ``10x10x10``:

.. code-block:: python

    nfield = nm.DirectorField(
        mesh_lengths=(10,10,10), mesh_dimensions=(80,80,80))

Next, we need to specify numerical values for the director field. Two methods are possible:
either you already have a numpy array containing the values of your director field, in which
case you can directly give this array to the ``DirectorField`` object (remember, you need to
make sure that this array is of shape ``(Nz,Ny,Nx,3)``):

.. code-block:: python

    nfield.vals = my_director_vals_numpy_array

Or you have an analytical formula for the director field, in which case you can define three
python functions and give these to the ``DirectorField`` object. In this tutorial, we will
assume the latter option and define the director field of a double twist cylinder:

.. code-block:: python

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
    nfield.init_from_funcs(nx,ny,nz)

If the analytical formula for the director components do not give a normalized director values,
you can still normalize manually the director values after importing them:

.. code-block:: python

    nfield.normalize()

Finally, you can apply geometric transformation to the director field with the methods
``rotate_90deg``, ``rotate_180deg`` and ``extend``,
as well as specify a non-trivial domain for the LC phase with the method ``set_mask``.
All these methods are documented in the API section of this wiki. Here, we will simply
demonstrate the capabilities of the ``DirectorField`` object by applying a 90° rotation around
the axis ``x``, extending the mesh in the ``xy`` plane with a scale factor of 2, and defining a
droplet mask centered on the mesh with a diameter equal to the mesh height:

.. code-block:: python

    nfield.rotate_90deg("x")
    nfield.extend(2,2)
    nfield.set_mask(mask_type="droplet")

Note that extending the mesh in the xy direction is essential if you define a non-trivial
LC mask, because you need to leave enough room for the optical fields to propagate around the LC
domain.

And that's it! If you want to save this director file to a XML VTK file (the standard format
used by the excellent visualisation software `Paraview <https://www.paraview.org/>`_), you can
add the following command to your script:

.. code-block:: python

    nfield.save_to_vti("double_twist_droplet")


You can import back the generated file in any script by directly constructing the DirectorField
object with the path to this file:

.. code-block:: python

    nfield = nm.DirectorField(vti_file="double_twist_droplet.vti")

This functionality is especially useful if generating the director field values takes a lot of
time.



.. _lcmat:

Defining a LCMaterial
---------------------

The next step is to define possible isotropic layers above the LC layer (which can distort the
optical fields on the focal plane), as well as the refractive indices of all the materials in the
sampe. Since our system here consists of a droplet embedded in another fluid, we need to specify
both extraordinay and ordinary indices for the LC droplet and the refractive index of the host
fluid. All these informations are stored in the class LCMaterial:

.. code-block:: python

    mat = nm.LCMaterial(
        director_field=nfield, ne=1.5, no=1.7, nhost=1.55)

Note that you can also specify refractive indices with a string expression depending on the
wavelength variable "lambda", in case you want to take into account the dispersivity of the
materials of you sample. 

We also want to add a glass plate above the sample and additional space for the host fluid
between the droplet and the glass plate:

.. code-block:: python

    mat.add_isotropic_layer(nlayer=1.55, thickness=5) # 5 µm space between the droplet and glass plate
    mat.add_isotropic_layer(nlayer=1.51, thickness=1000) # 1mm-thick glass plate

We don't specify isotropic layers below the sample because the high-level interface only support
input optical fields propagating in the ``z`` direction (in which case the amplitude of the
fields is uniformly affected by any isotropic layers orthogonal to ``z``). This may change in
the future, since the backend ``dtmm`` does support multiple plane-wave source as in a real
Köhler illumination setup.

.. _prop:

Propagating optical fields through the sample
---------------------------------------------

Now that the sample geometry is fully caracterized, we can propagate fields through the sample
and back to the central focal plane. This is simple as defining an array of wavelengths defining
the spectrum extent of the light source, creating a ``LightPropagator`` object, and calling the 
method ``propagate_field``:

.. code-block:: python

    wavelengths = np.linspace(0.4, 0.8, 11)
    sim = nm.LightPropagator(
        material=mat, wavelengths=wavelengths, numerical_aperture=0.4)
    output_fields = sim.propagate_fields(method="bpm")

The numerical aperture defined in this code snippet corresponds to the one of the microscope
objective. The ``propagate_fields`` method uses the specified backend to propagate fields (here,
``bpm-solver``) and returns an ``OpticalFields`` object containing the results of the simulation.
Periodic boundary conditions in the ``x`` and ``y`` directions are always assumed, so you should
always extend apropriately your director field in order to have a uniform field near the mesh
boundaries.

Note that internally two simulations are run, one with an input light source polarised along
``x`` and the other with an input light source polarised along ``y``. This allows us to
fully caracterize the transmission matrix of the sample and reconstruct any type of micrographs
(bright field, crossed polariser...), as we will see in the next section.

Similaryly to the ``DirectorField`` object, you can save the output fields to a XML VTK file, and
reimport them in other scripts:

.. code-block:: python

    # If you want to save the simulation results
    output_fields.save_to_vti("optical_fields")

    # If you want to reimport saved simulation results
    output_fields = nm.OpticalFields(vti_file="optical_fields.vti")


.. _viz:

Visualising optical micrographs
-------------------------------

To help the user visualise optical micrographs as in a real microscope, ``nemaktis`` includes
a graphical user interface allowing to generate any type of micrograph in real-time. Once
you have generated/imported optical fields in you script, you can start using this interface
with the following lines of code:

.. code-block:: python

    viewer = nm.FieldViewer(output_fields)
    viewer.plot()

All parameters in this user interface should be pretty self-explanatory. We will simply mention
here that the quarter and half wavelengths compensators are assumed to be achromatic, while the
full-wave "tint sensitive" compensator is aproximated with a slab of wavelength-independent
refractive index with a full-wave shift only at a wavelength of 540 nm.

Concerning color management, we assume a D65 light source and project the output light spectrum
first on the XYZ space, then on the sRGB color space, to finally obtain a usual RGB picture. 
For more details, see `<https://dtmm.readthedocs.io/en/latest/tutorial.html#color-conversion>`_.

Finally, refocalisation of the optical micrographs is done by switching to Fourrier space and
using the exact propagator for the Helmholtz equation in free space.
