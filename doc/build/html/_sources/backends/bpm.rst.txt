Beam propagation backend
========================

The executable of the beam propagation backend is named ``bpm-solver`` and
can be called in any terminal with the conda environment for nemaktis
activated.

``bpm-solver`` relies on json-like configuration file, so you don't have
to touch to any C++ code if you just want to use the method. To generate
a default configuration file, run:

.. code-block::
    
    bpm-solver -c [Name of the configuration file]

All parameters in this configuration file are fully documented, so you
should be able to understand how to make the code working just by
reading and ajusting the parameters in this file. More information on
the subtleties of this code will be added later on this wiki, for now we
will just mention that the input vti file for the director field can be
created directly with the high-level interface (see
:class:`~nemaktis.lc_material.DirectorField`).

Once you are satisfied with the change you made to your configuration
file, you can actually run the code by typing:

.. code-block::

    bpm-solver -x [Name of the configuration file]

