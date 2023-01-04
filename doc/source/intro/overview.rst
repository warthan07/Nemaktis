.. _overview:

Overview
========

Nemaktis is an open-source platform including tools for propagating and visualising optical
fields in complex birefringent media such as liquid crystal (LC) layers. It includes three
backends implementing advanced numerical methods for light propagation, as well as an
easy-to-use high level interface in python allowing to quickly setup a simulation and visualize
optical microraphs of a LC structure as in a real microscope. It goes well beyond the
Jones method usually used in LC research, by accurately modeling diffraction, walk-off,
focusing effects, Koehler illumination... 

If you want a platform for easily comparing experimental images and numerical micrographs from
simulated or theoretical birefringent structures, you are in the right place!


License
-------

Nemaktis is released under MIT license so you can use it freely. Please cite the following
publications if you use Nemaktis to prepare a figure in a scientific paper:

* `G. Poy, S. Žumer, Soft Matter 15, 3659-3670 (2019) <https://doi.org/10.1039/C8SM02448K>`_
* `G. Poy, S. Žumer, Optics Express 16, 24327-24342 (2020) <https://doi.org/10.1364/OE.400984>`_


Contributors
------------

* High-level interface, ray-tracing and beam propagation backends: Guilhem Poy, Slobodan Žumer.
* Diffraction transfer matrix backend: Andrej Petelin, Alex Vasile.


Highlights
----------

* Easy-to-use scripting interface in python
* Support for Koehler illumination setup (multiple incoming plane waves)
* Support for arbitrary number of isotropic layers around the birefringent object (e.g.
  glass plates).
* Support for arbitrary uniaxial media (biaxial support coming soon).
* Graphical user interface to visualize optical fields, with interactive sliders for the
  parameters of the microscope.


Limitations
-----------

* Paraxial propagation is assumed for the BPM backend (but the PSF of the microscope can be
  set with a high NA).
* No support for reflection microscopy
* No support (yet) for biaxial media

The last two limitations should be lifted in the future since the associated theoretical
framework is already ready (I just need to find students or the time to code everything!).
If you want to implement new features that you think could benefit the whole software,
please contact me (address below)!

Concerning the paraxial approximation, I can mention that I have also developped a
closed-source wide-angle beam propagation method, which can model non-linear optical
systems, wide-angle deflection of light beams by birefringent structures, waveguiding...
If you are interested in starting a collaboration on this closed-source software, please
send me a quick message explaining the optics problem that you want to solve.


Contact
-------

* Personal website: https://www.syncpoint.fr
* Email: guilhempoy07 [*at*] gmail [*dot*] com
