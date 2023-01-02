__version__ = "1.4.6"

import matplotlib
matplotlib.use("Qt5Agg")

from .lc_material import LCMaterial, TensorField, QTensorField, DirectorField
from .light_propagator import LightPropagator, OpticalFields
from .field_viewer import FieldViewer
