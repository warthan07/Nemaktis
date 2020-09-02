.. _microscopy_model:

Microscopy model for Nemaktis
=============================

We present here the theoretical model of microscopy that is at the core of Nemaktis. A few
interactive graphics are provided in order to better understand important concepts. The
javascript code for these interactive examples can be found in an `ObservableHQ notebook
<https://observablehq.com/@warthan07/microscopy-model-for-nemaktis>`_.

General description
-------------------

In a real transmission and/or reflection microscope, objects are imaged using a combination
of lighting systems and lenses. The path of light in such microscopes can always be
decomposed in three steps:

1. Light propagates from the illuminating source to the object through the illumination
   optical setup.

2. Light is transmitted through (or reflected from) the object.

3. Light transmitted or reflected from the object propagates through the microscope
   objective and form an observable image in a target imaging plane.

The case of spectrally-extended lighting (e.g. white light lamp) can be easily covered by
projecting on an appropriate color space the final images formed by the different
wavelengths of the lamp spectrum. In Nemaktis, this is done internally after repeating the
imaging simulations for all the wavelengths in a user-defined array approximating the lamp
spectrum. For more details on the color space projection method, see `Color
conversion <https://dtmm.readthedocs.io/en/latest/tutorial.html#color-conversion>`_ in the
documentation of ``dtmm``, one of the backend used in Nemaktis.  Here, we consider for
simplicity's sake a simple microscopy model based on lighting with a single wavelength. We
describe in the next sections the physical mechanisms behind the three steps introduced
above, as schematized below in a simplified representation of our virtual microscope in
transmission mode:

.. raw:: html

  <div id="microscope-fig">
    <div class="observablehq-chart_microscope"></div>
  </div>
  <script type="module">
    import {getRuntime} from "../_static/observable.js"
    import {Inspector} from "https://cdn.jsdelivr.net/npm/@observablehq/runtime@4/dist/runtime.js";
    import notebook from "https://api.observablehq.com/@warthan07/microscopy-model-for-nemaktis.js?v=3";
    getRuntime("#microscope-fig").module(notebook, name => {
      if(name === "chart_microscope") return Inspector.into("#microscope-fig .observablehq-chart_microscope")();
    });
  </script>


Koehler illumination setup
--------------------------

The first propagation step is the illumination of the object by the light source. The
standard illumination setup used in most microscopes is called the Koehler illumination
setup (introduced by August Koehler in 1893), and has the advantage of allowing a uniform
lighting even with non-uniform light source. In short, it allows to map each point of the
light source to a single uniform plane wave incident on the object with a certain angle; the
maximum angle of incidence for the plane waves is set by an aperture called the **condenser
aperture**, thus the set of plane waves incident on the object all have wavevectors included
inside a cone of illumination whose opening is set by the condenser aperture. In addition, a
**field aperture** allows to control the size of the lighting spot on the object, but this
feature is irrelevant in Nemaktis since we always assume that the whole computational box
representing the object is illuminated.

In order to better understand how this illumination setup works, an interactive example is
provided below, where the reader can dynamically adjust the sliders for opening/closing the
field and condenser apertures:

.. raw:: html

  <div id="koehler-fig">
    <div class="observablehq-viewof-cond_ap_opening"></div>
    <div class="observablehq-viewof-field_ap_opening"></div>
    <div class="observablehq-chart_koehler"></div>
  </div>
  <script type="module">
    import {getRuntime} from "../_static/observable.js"
    import {Inspector} from "https://cdn.jsdelivr.net/npm/@observablehq/runtime@4/dist/runtime.js";
    import notebook from "https://api.observablehq.com/@warthan07/microscopy-model-for-nemaktis.js?v=3";
    getRuntime("#koehler-fig").module(notebook, name => {
      if(name === "viewof cond_ap_opening") return Inspector.into("#koehler-fig .observablehq-viewof-cond_ap_opening")();
      if(name === "viewof field_ap_opening") return Inspector.into("#koehler-fig .observablehq-viewof-field_ap_opening")();
      if(name === "chart_koehler") return Inspector.into("#koehler-fig .observablehq-chart_koehler")();
    });
  </script>

A correctly assembled Koehler illumination setup has the following properties:

* The field aperture is in the back focal plane of the lamp collector lens.
* The condenser aperture is in the front focal plane of the condenser lens.
* The image of the lamp filament through the lamp collector lens is in the same plane as the
  condenser aperture.
* The image of the field aperture throught the condenser lens is is in the same plane as the
  object.

We emphasize that the lamp filament is always spatially incoherent, thus the different
incident plane waves cannot interfer between themselves. This means that the final image in
the microscope is always obtained by summing-by-intensity the individual images formed by
each incident plane waves. In real life, there is always an infinite number of plane waves
incident on the object, but in the computer one must choose an approximate discrete set of
plane waves. In Nemaktis, the set of incoming plane waves is chosen to have the following
wavevectors (assuming that the third coordinate correspond to the main propagation axis in
the microscope):

.. math::

  \vec{k}^{(k,l)}=k_0\left(\begin{aligned}
    q^{(k)} \cos\theta^{(k,l)} \\ q^{(k)} \sin\theta^{(k,l)} \\ \sqrt{1-\left[q^{(k)}\right]^2}
  \end{aligned}\right)

where we defined :math:`k_0=2\pi/\lambda` with :math:`\lambda` the wavelength in empty space and:

.. math::

  \begin{aligned}
    q^{(k)} &= \frac{k}{N_r-1}\mathrm{NA}_\mathrm{max},\quad\quad k=0\cdot\cdot\cdot N_r-1 \\
    \theta^{(k,l)} &= \frac{\pi l}{3k},\quad\quad\quad\quad\quad\quad l=0\cdot\cdot\cdot 6k
  \end{aligned}

Here, :math:`\mathrm{NA}_\mathrm{max}=\sin\psi_\mathrm{max}` (with :math:`\psi_\mathrm{max}`
the maximal angle of opening of the wavevectors) is the maximal numerical aperture of the
Koehler illumination setup, and :math:`N_r` correspond to the number of discretization steps
in the radial direction. This choice of wavevectors correspond to a standard discretization
of a circular aperture in the transverse plane, which can be interactively visualized by
adjusting the sliders for :math:`N_r` and :math:`\mathrm{NA}`.
