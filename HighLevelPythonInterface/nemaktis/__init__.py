__version__ = "2.0"

import matplotlib
matplotlib.use("Qt5Agg")

from .lc_material import TensorField, DirectorField, QTensorField, MicroscopeSample
from .light_propagator import LightPropagator, OpticalFields
from .field_viewer import FieldViewer
