# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#

from unittest.mock import MagicMock, Mock, patch
import sys, os, six
sys.path.insert(0, os.path.abspath("../../HighLevelPythonInterface/"))

# Modules that we need to mock
MOCK_MODULES = [
    "traits",
    "traits.api",
    "traits.etsconfig",
    "traits.etsconfig.api",
    "traitsui",
    "traitsui.api",
    "traitsui.qt4",
    "traitsui.qt4.editor",
    "traitsui.qt4.basic_editor_factory",
    "traitsui.delegating_handler",
]

# Collect the types from traitsui that are based on HasTraits
# We will need to mock them in a special way so that
# TraitDocumenter can properly identify and document traits.
#  from traits.api import HasTraits, HasPrivateTraits

class HasTraits(MagicMock):
    @classmethod
    def __getattr__(cls, name):
        return MagicMock()

class HasPrivateTraits(MagicMock):
    @classmethod
    def __getattr__(cls, name):
        return MagicMock()

MOCK_TYPES = []
MOCK_TYPES.append(
    ("traitsui.delegating_handler", "DelegatingHandler", (HasTraits,))
)
MOCK_TYPES.append(
    ("traitsui.view_element", "ViewSubElement", (HasPrivateTraits,))
)
MOCK_TYPES.append(
    ("traitsui.editor", "Editor", (HasPrivateTraits,))
)
MOCK_TYPES.append(
    ("traitsui.basic_editor_factory", "BasicEditorFactory", (HasPrivateTraits,))
)


# Create the custom types for the HasTraits based traitsui objects.
TYPES = {
    mock_type: type(mock_type, bases, {"__module__": path})
    for path, mock_type, bases in MOCK_TYPES
}

class DocMock(MagicMock):
    @classmethod
    def __getattr__(self, name):
        if name in ("__file__", "__path__", "_mock_methods"):
            # sphinx does not like getting a Mock object in this case.
            return "/dev/null"
        else:
            # Return a mock or a custom type as requested.
            return TYPES.get(name, DocMock(mocked_name=name))

    # MagicMock does not define __call__ we do just to make sure
    # that we cover all cases.
    def __call__(self, *args, **kwards):
        return DocMock()

    @property
    def __name__(self):
        # Make sure that if sphinx asks for the name of a Mocked class
        # it gets a nice strings to use (instead of "DocMock")
        return self.mocked_name

# Add the mocked modules to sys
sys.modules.update(
    (mod_name, DocMock(mocked_name=mod_name)) for mod_name in MOCK_MODULES
)


# -- Project information -----------------------------------------------------

project = 'Nemaktis'
copyright = '2019, Guilhem Poy, Andrej Petelin'
author = 'Guilhem Poy, Andrej Petelin'

# The full version, including alpha/beta/rc tags
release = '1.0.2'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
pygments_style = 'sphinx'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []

autodoc_mock_imports = ['pyface', 'dtmm', 'matplotlib', 'scipy', 'numpy',
        'vtk', 'bpm_backend', 'pyfftw']
autodoc_member_order = 'bysource'
    
