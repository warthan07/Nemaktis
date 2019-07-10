__version__ = "0.1.0.dev"

import matplotlib
matplotlib.use("Qt5Agg")

from .lc_material import LCMaterial, DirectorField
from .light_propagator import LightPropagator, OpticalFields
from .field_viewer import FieldViewer
