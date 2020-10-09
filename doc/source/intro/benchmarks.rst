.. _benchmarks:

Accuracy and efficiency of Nemaktis' backends
=============================================

We present here accuracy and efficiency benchmarks for the main two backends of Nemaktis:
DTMM and BPM. The javascript code for the interactive graphs of this page can be found in an
`ObservableHQ notebook
<https://observablehq.com/@warthan07/accuracy-and-efficiency-of-nemaktis-backends>`_

1. Introduction
---------------

As explained in `[Transmission/Reflection of light inside the object]
<https://nemaktis.readthedocs.io/en/latest/intro/microscopy_model.html#transmission-reflection-of-light-inside-the-object>`_,
the two main backends of Nemaktis (DTMM and BPM) propagate optical fields through a series
of inhomonogeneous birefringent layers based on different formulation of the evolution
operator of Maxwell equations. Both schemes include diffraction effects due to
inhomogeneities in the optical and/or permittivity fields. Since analytical solutions of
Maxwell equations are readily available only on homogeneous birefringent media, we propose
here to evaluate the accuracy of both schemes by simulating the diffraction of a
highly-focused Gaussian beam in a uniform uniaxial crystal. This has the advantage of
evaluating in one simulation the propagation accuracy of many different plane waves since a
focused Gaussian beam spans a wide-area in the transverse Fourier plane.

However, the reader should keep in mind that "real-life" computational errors of simulations
on inhomogeneous birefringent media can be higher than the ones presented here, since they
can also include contributions from permittivity-based diffraction (e.g. the boundary of a
birefringent droplet, or sharp variations of the optical axis). In general, one should
always try to adapt the mesh resolution in order to resolve the finest details of the
permittivity tensor field and get the best possible accuracy. Note that discontinuities of
permittivity are accurately modeled in Nemaktis if they are positioned midway between two
mesh points, which is why fine meshes are needed to model complex boundaries.

Our methodology to calculate the computational error is as follows:

1. Define an input optical field with polarisation :math:`\vec{u}`:

.. math::

  \vec{E}(\vec{r}_\perp,z=0)=\exp\left[-\frac{\left|\vec{r}_\perp\right|^2}{2w^2}\right]\vec{u}

2. Propagate the optical field on a computational box with length L for each side.

3. Calculate the mean error associated with all Fourier components relevant for microscopy
   simulations (i.e. paraxial components):

.. math::

  \epsilon=\left\langle\frac{\left|
  \tilde{\vec{E}}\left[\vec{k}_\perp,z\right]-\tilde{\vec{E}}_{\rm exact}\left[\vec{k}_\perp,z\right]
  \right|}{2\left|\tilde{\vec{E}}_{\rm exact}\left[\vec{k}_\perp,z\right]\right|}
  \right\rangle_{\theta(\vec{k}_\perp)<\theta_m}

In the definition above, a tilde indicate a partial Fourier transform along the :math:`x`
and :math:`y` coordinates (associated with a Fourier frequency :math:`\vec{k}_\perp`),
:math:`\theta(\vec{k}_\perp)` is the the angle between :math:`\vec{k}_\perp` and the main
direction of propagation :math:`z`, and the bracket indicate a constrained average on
paraxial Fourier frequencies and :math:`z` coordinates. The exact solution
:math:`\tilde{\vec{E}}_{\rm exact}\left[\vec{k}_\perp,z\right]` is obtained directly from an
eigenmode decomposition of Maxwell equations in transverse Fourier space, which can be done
analytically since we assumed a uniform optical axis (see python script `propagate_fields.py
<https://github.com/warthan07/Nemaktis/blob/master/HighLevelPythonInterface/examples/benchmarks/3D/propagate_fields.py>`).
The maximum angle of propagation :math:`\theta_m\approx24°` was choosen based on the typical
numerical aperture :math:`\sin\theta_m=0.4` of microscope objectives.

Note that this definition of the error is roughly independent from the choice of input
profile, since it basically measure the relative error made on phase evolution in Fourier
space (hence the factor 2: a phase error of :math:`\pi` is associated with a maximal error
of 1). However, it still necessitates a focused profile in order to get non-negligible
high-frequency components and avoid division by zero in the definition of :math:`\epsilon`. 

In addition to the calculation of the error, we also sytematically save the computational
time associated with one calculation. We chose a waist :math:`w=2` µm and a mesh length
:math:`L=10` µm in all the numerical experiments of this page, and vary the mesh density and
choice of backends (BPM or DTMM(D) schemes, where D represents the "diffraction" parameter).
You can experimentally test the accuracy of our schemes using the scripts provided in
the `benchmarks folder of the GitHub repository
<https://github.com/warthan07/Nemaktis/tree/master/HighLevelPythonInterface/examples/benchmarks>`_

2. Results with 2D meshes
-------------------------

We present first our results on 2D meshes (i.e. mesh with one transverse direction and one
longitudinal direction z associated with the main direction of propagation of the Gaussian
beam). Below, we plot the computational error as a function of the mesh spacing for
different choice of backends:`

.. raw:: html

  <div id="err-2D-fig">
    <div class="observablehq-err_vs_d_chart_2D"></div>
  </div>
  <script type="module">
    import {getRuntime} from "../_static/observable.js"
    import {Inspector} from "https://cdn.jsdelivr.net/npm/@observablehq/runtime@4/dist/runtime.js";
    import notebook from "https://api.observablehq.com/@warthan07/accuracy-and-efficiency-of-nemaktis-backends.js?v=3";
    getRuntime("#err-2D-fig").module(notebook, name => {
      if(name === "err_vs_d_chart_2D") return Inspector.into("#err-2D-fig .observablehq-err_vs_d_chart_2D")();
    });
  </script>

The saturation of the computational mesh on very fine mesh can be simply understood from the
fact that all backends in Nemaktis are approximate and not exact, even with an infinite
number of degree of freedom. However, they still allows to get very good computational
error, which can be below 1% for the BPM, DTMM(5) and DTMM(7) backends. The inacurracy of
the DTMM(1) and DTMM(3) backends is due to the fact that these backends are unable to
correctly propagate tilted plane waves (i.e. high frequency transverse Fourier modes), so if
diffraction is not negligible in the system you want to study (tight beam profile or
birefringent structures with details of typical length 1-10 µm) you should avoid these
schemes.

As can be observed, the BPM backend allows to get the best possible accuracy provided a
sufficiently fine mesh is used. However, the DTMM schemes have a very interesting property:
they are less sensitive to the mesh spacing than the BPM scheme. This means that you will
always get a smaller computational error on small meshes by using high-order DTMM schemes
instead of the BPM accuracy.

But accuracy is not the only interesting parameter here, since the running times of these
simulations widely vary when changing the backend and/or mesh spacing. To understand how all
these parameters are linked, we provide below an interactive bar chart allowing to quickly
visualize the error and running times of each backend for a given mesh spacing. Since 2D
simulations are computationally not very intensive, the running times were evaluated on a 6
years old laptop with a processor i7-4600M (4 threads).

.. raw:: html

  <div id="times-2D-fig">
    <div class="observablehq-viewof-dy_idx_2D"></div>
    <div class="observablehq-viewof-order_by_2D"></div>
    <div class="observablehq-err_times_chart_2D"></div>
    <div class="observablehq-err_times_chart_2D_update"></div>
  </div>
  <script type="module">
    import {getRuntime} from "../_static/observable.js"
    import {Inspector} from "https://cdn.jsdelivr.net/npm/@observablehq/runtime@4/dist/runtime.js";
    import notebook from "https://api.observablehq.com/@warthan07/accuracy-and-efficiency-of-nemaktis-backends.js?v=3";
    getRuntime("#times-2D-fig").module(notebook, name => {
      if(name === "viewof dy_idx_2D") return Inspector.into("#times-2D-fig .observablehq-viewof-dy_idx_2D")();
      if(name === "viewof order_by_2D") return Inspector.into("#times-2D-fig .observablehq-viewof-order_by_2D")();
      if(name === "err_times_chart_2D") return Inspector.into("#times-2D-fig .observablehq-err_times_chart_2D")();
      if(name === "err_times_chart_2D_update") return Inspector.into("#times-2D-fig .observablehq-err_times_chart_2D_update")();
    });
  </script>

Not very surprinsingly, the inacurate DTMM(1) and DTMM(3) backends are also the fastest.
Basically, these low-order DTMM schemes correspond to Jones-like calculus with a
fast-but-inacurate treatment of diffraction, which is why their computational error is high
due to the presence of high-frequency Fourier modes in these simulations. But if you know in
advance that diffraction in your system is negligible (for example if the optical axis vary
over lengths much bigger than the wavelength along directions orthogonal to the main axis of
propagation), these schemes are a really good choice since they are very fast and can still
be reasonably accurate for propagating low-frequency Fourier modes.

As for the BPM, DTMM(3) and DTMM(5), it can be observed that the DTMM schemes wins the time
race on small meshes, while the BPM schemes is the fastest (and most accurate) on big
meshes. This can be interpreted from the complexity of the numerical algorithm of each
backends: the DTMM(D) backend has a :math:`O\left(D^{(d-1)} N \log\left[N/N_z\right]\right)`
complexity while the BPM backend has a better linear complexity :math:`O(N)`, with :math:`d`
the dimensionality of the mesh, :math:`N` the total number of mesh point and :math:`N_z` the
number of points along the z axis); therefore, it is not really surprising that the DTMM
schemes gets penalized in terms of running times for high diffraction parameter :math:`D` or
points number :math:`N`.


3. Results with 3D meshes
-------------------------

We now turns our focus to 3D meshes (i.e. meshes with two transverse directions and one
longitudinal direction z associated with the main direction of propagation of the Gaussian
beam). Results are qualitatively similar than for 2D meshes, except now the DTMM(7) backend
is always the most accurate scheme, whatever the mesh spacing. Nevertheless, the DTMM(5) and
BPM backends still manage to get relatively good computational error of :math:`\sim` 1% on
sufficiently fine mesh.`

.. raw:: html

  <div id="err-3D-fig">
    <div class="observablehq-err_vs_d_chart_3D"></div>
  </div>
  <script type="module">
    import {getRuntime} from "../_static/observable.js"
    import {Inspector} from "https://cdn.jsdelivr.net/npm/@observablehq/runtime@4/dist/runtime.js";
    import notebook from "https://api.observablehq.com/@warthan07/accuracy-and-efficiency-of-nemaktis-backends.js?v=3";
    getRuntime("#err-3D-fig").module(notebook, name => {
      if(name === "err_vs_d_chart_3D") return Inspector.into("#err-3D-fig .observablehq-err_vs_d_chart_3D")();
    });
  </script>

However the running times of DTMM backends vs BPM backend are vastly different than in the
2D case, as expected from the :math:`D^{(d-1)}` factor in the complexity of DTMM backends
(see above). Since 3D simulations are computationally intensive, the results of the
interactive bar chart below were obtained on a recent desktop computer with processor
i7-7800X (12 threads).

.. raw:: html

  <div id="times-3D-fig">
    <div class="observablehq-viewof-dy_idx_3D"></div>
    <div class="observablehq-viewof-order_by_3D"></div>
    <div class="observablehq-err_times_chart_3D"></div>
    <div class="observablehq-err_times_chart_3D_update"></div>
  </div>
  <script type="module">
    import {getRuntime} from "../_static/observable.js"
    import {Inspector} from "https://cdn.jsdelivr.net/npm/@observablehq/runtime@4/dist/runtime.js";
    import notebook from "https://api.observablehq.com/@warthan07/accuracy-and-efficiency-of-nemaktis-backends.js?v=3";
    getRuntime("#times-3D-fig").module(notebook, name => {
      if(name === "viewof dy_idx_3D") return Inspector.into("#times-3D-fig .observablehq-viewof-dy_idx_3D")();
      if(name === "viewof order_by_3D") return Inspector.into("#times-3D-fig .observablehq-viewof-order_by_3D")();
      if(name === "err_times_chart_3D") return Inspector.into("#times-3D-fig .observablehq-err_times_chart_3D")();
      if(name === "err_times_chart_3D_update") return Inspector.into("#times-3D-fig .observablehq-err_times_chart_3D_update")();
    });
  </script>

This time, the BPM backend is practically always faster than DTMM schemes (only the DTMM(1)
can be faster than BPM on fine meshes), while having a very good computational error for
most mesh spacings. In particular, the very accurate DTMM(5) and DTMM(7) schemes
necessitates 4-8 times longer running times than the BPM scheme. As a consequence, we
recommend to use the DTMM schemes on 3D meshes only when you want a fast simulation method
without accurate diffraction (DTMM(1) backend) or a very accurate but very slow simulation
(DTMM(7) backend). For all other case of applications, the BPM backend provide a reliable
and accurate simulation scheme whatever the size of the computational mesh.
