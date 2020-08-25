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
   objective and form an observable image in a target plane.

The case of spectrally-extended lighting (e.g. white light lamp) can be easily covered by
projecting on an appropriate color space the final images formed by the different
wavelengths of the lamp spectrum. In Nemaktis, this is done internally after repeating the
imaging simulations for all the wavelengths in a user-defined array approximating the lamp
spectrum. For more details on the color space projection method, see `Color
conversion<https://dtmm.readthedocs.io/en/latest/tutorial.html#color-conversion>`_ in the
documentation of ``dtmm``, one of the backend used in Nemaktis.  Here, we consider for
simplicity's sake a simple microscopy model based on lighting with a single wavelength. We
describe below the physical mechanisms behind the three steps introduced above.


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
    import {getRuntime} from "_static/observable.js"
    import {Inspector} from "https://cdn.jsdelivr.net/npm/@observablehq/runtime@4/dist/runtime.js";
    import notebook from "https://api.observablehq.com/@warthan07/microscopy-model-for-nemaktis.js?v=3";
    getRuntime("koehler-fig").module(notebook, name => {
      if(name === "viewof cond_ap_opening") return Inspector.into("#koehler-fig .observablehq-viewof-cond_ap_opening")();
      if(name === "viewof field_ap_opening") return Inspector.into("#koehler-fig .observablehq-viewof-field_ap_opening")();
      if(name === "chart_koehler") return Inspector.into("#koehler-fig .observablehq-chart_koehler")();
    });
  </script>
