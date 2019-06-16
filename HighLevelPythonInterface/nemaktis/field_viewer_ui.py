import six
from traits.etsconfig.api import ETSConfig
ETSConfig.toolkit = 'qt4'

import matplotlib
matplotlib.use("Qt4Agg", warn=False, force=True)
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvas

from traits.api import *
from traitsui.api import *
from traitsui.qt4.editor import Editor
from traitsui.qt4.basic_editor_factory import BasicEditorFactory

from pyface.qt import QtGui, QtCore
from pyface.util.guisupport import start_event_loop_qt4

from threading import Thread
from time import sleep
import os


class OpticalElementSettings(HasTraits):
    polariser = Enum(("Yes", "No"))
    analyser = Enum(("Yes", "No"))
    compensator = Enum(("No", "Quarter-wave", "Half-wave", "Tint-sensitive"))

    polariser_angle = Float(0)
    analyser_angle = Float(90)
    compensator_angle = Float(0)

    view = View(
        Group(
            Spring(),
            Group(
                Item("polariser", label="Polariser?", style="custom"),
                Item("compensator", label="Compensator?"),
                Item("analyser", label="Analyser?", style="custom"),
                Spring(width=-50),
                Spring(width=-50),
                Spring(width=-50),
                Item("polariser_angle", enabled_when="polariser==\"Yes\"", editor=RangeEditor(
                    mode="slider", format="%.2f",
                    low=-90., high=90., low_label="-90", high_label="90")),
                Item("compensator_angle", enabled_when="compensator!=\"No\"", editor=RangeEditor(
                    mode="slider", format="%.2f",
                    low=-90., high=90., low_label="-90", high_label="90")),
                Item("analyser_angle", enabled_when="analyser==\"Yes\"", editor=RangeEditor(
                    mode="slider", format="%.2f",
                    low=-90., high=90., low_label="-90", high_label="90")),
                orientation="horizontal",
                columns=3),
            Spring(),
            orientation="vertical",
            label="Optical elements settings",
            show_border = True),
        resizable=True)

    def __init__(self, callback_dict):
        self._callback_dict = callback_dict

    def _polariser_changed(self):
        self._callback_dict["set_polariser"](self.polariser)
    def _analyser_changed(self):
        self._callback_dict["set_analyser"](self.analyser)
    def _compensator_changed(self):
        self._callback_dict["set_compensator"](self.compensator)
    def _polariser_angle_changed(self):
        self._callback_dict["set_polariser_angle"](self.polariser_angle)
    def _analyser_angle_changed(self):
        self._callback_dict["set_analyser_angle"](self.analyser_angle)
    def _compensator_angle_changed(self):
        self._callback_dict["set_compensator_angle"](self.compensator_angle)



class MicroscopeSettings(HasTraits):
    intensity = Float(1)
    min_intensity = Float(0)
    max_intensity = Float(2)

    z_focus = Float(0)
    min_focus = Float(-10)
    max_focus = Float(10)

    color_mode = Enum(("RGB", "Grayscale"))
    n_tiles_x = Int(1)
    n_tiles_y = Int(1)

    view = View(
        Group(
            Spring(),
            Group(
                Item("intensity", editor=RangeEditor(
                    mode="slider", format="%.2f", is_float=True,
                    low_name="min_intensity", high_name="max_intensity")),
                Item("z_focus", label="Z-focus", editor=RangeEditor(
                    mode="slider", format="%.2f", is_float=True,
                    low_name="min_focus", high_name="max_focus")),
                Spring(width=-50),
                Spring(width=-50),
                Item("min_intensity", label="Min", width=-40),
                Item("min_focus", label="Min", width=-40),
                Spring(width=-50),
                Spring(width=-50),
                Item("max_intensity", label="Max", width=-40),
                Item("max_focus", label="Max", width=-40),
                orientation="horizontal",
                columns=2),
            Spring(),
            Group(
                Item("n_tiles_x", width=-40, editor=RangeEditor(
                    mode="spinner", low=1, high=10)),
                Spring(width=-30),
                Item("n_tiles_y", width=-40, editor=RangeEditor(
                    mode="spinner", low=1, high=10)),
                Spring(),
                orientation="horizontal"),
            Spring(),
            Group(
                Item("color_mode", label="Color mode", editor=EnumEditor(
                    values={'RGB': 'RGB', 'Grayscale': 'Grayscale'},
                    format_func=six.text_type, cols=2), style="custom")),
            Spring(),
            orientation="vertical",
            show_border=True,
            label="Microscope settings"),
        resizable=True)

    def __init__(self, callback_dict):
        self._callback_dict = callback_dict

    def _intensity_changed(self):
        self._callback_dict["set_intensity"](self.intensity)
    def _z_focus_changed(self):
        self._callback_dict["set_z_focus"](self.z_focus)
    def _color_mode_changed(self):
        self._callback_dict["set_color_mode"](self.color_mode)
    def _n_tiles_x_changed(self):
        self._callback_dict["set_n_tiles_x"](self.n_tiles_x)
    def _n_tiles_y_changed(self):
        self._callback_dict["set_n_tiles_y"](self.n_tiles_y)


class TextDisplay(HasTraits):
    string =  String()
    view= View(Item("string", show_label=False, springy=True, style="readonly"))


class TempDisplayThread(Thread):
    def __init__(self, display, message, duration):
        super(TempDisplayThread, self).__init__()
        self.wants_abort = False
        self.display = display
        self.message = message
        self.duration = duration

    def run(self):
        # Blink effect in case there is already something on the display
        self.display.string = ""
        sleep(0.1)

        self.display.string = self.message
        time = 0.
        while not self.wants_abort and time<self.duration:
            sleep(0.1)
            time += 0.1
        self.display.string = ""


class SaveSettings(HasTraits):
    directory = Directory(os.getcwd())
    base_name = Str("micrograph")
    save = Button("Save to png")
    display = Instance(TextDisplay, ())
    _temp_display_thread = Instance(TempDisplayThread)

    view = View(
        Group(
            Spring(),
            Item("directory"),
            Item("base_name"),
            Group(
                Item("save", width=-100, show_label=False),
                Item("display", style="custom", show_label=False),
                orientation="horizontal"),
            Spring(),
            show_border=True,
            label="Save current micrograph"),
        resizable=True)

    def __init__(self, callback_dict):
        self._save_callback = callback_dict["save"]

    def _save_fired(self):
        try:
            self._save_callback(self.directory+"/"+self.base_name+".png")
        except Exception as e:
            if self._temp_display_thread and self._temp_display_thread.isAlive():
                self._temp_display_thread.wants_abort = True
                self._temp_display_thread.join()
            self._temp_display_thread = TempDisplayThread(self.display, str(e), 5)
            self._temp_display_thread.start()
        else:
            if self._temp_display_thread and self._temp_display_thread.isAlive():
                self._temp_display_thread.wants_abort = True
                self._temp_display_thread.join()
            self._temp_display_thread = TempDisplayThread(self.display, "Done!", 2)
            self._temp_display_thread.start()


class SettingPanels(HasTraits):
    optical_elements_settings = Instance(OpticalElementSettings)
    microscope_settings = Instance(MicroscopeSettings)
    save_settings = Instance(SaveSettings)

    view = View(
        Group(
            Spring(),
            Item("optical_elements_settings",
                style='custom', show_label=False, height=-140),
            Spring(height=-20),
            Item("microscope_settings",
                style='custom', show_label=False, height=-200),
            Spring(height=-20),
            Item("save_settings",
                style='custom', show_label=False, height=-140),
            Spring(),
            orientation="vertical"),
        resizable=True)

    def __init__(self, callback_dict):
        self.optical_elements_settings = OpticalElementSettings(callback_dict)
        self.microscope_settings = MicroscopeSettings(callback_dict)
        self.save_settings = SaveSettings(callback_dict)


class MPLFigureEditor(Editor):
    scrollable  = True

    def init(self, parent):
        self.control = self._create_canvas(parent)
        self.set_tooltip()

    def update_editor(self):
        pass

    def _create_canvas(self, parent):
        return FigureCanvas(self.value)


class MPLFigureEditorFactory(BasicEditorFactory):
    klass = MPLFigureEditor


class FigurePanel(HasTraits):
    figure = Instance(Figure, ())

    view = View(
        Group(
            Spring(height=-20),
            Item("figure", editor=MPLFigureEditorFactory(), show_label=False),
            Spring(height=-20),
            orientation="vertical"),
        resizable=True)

    def __init__(self, figure):
        super(FigurePanel, self).__init__()
        self.figure = figure


class FieldViewerUI(HasTraits):
    setting_panels = Instance(SettingPanels)
    figure_panel = Instance(FigurePanel)

    view = View(
        Group(
            Spring(width=-20),
            Item("figure_panel", style="custom", show_label=False),
            Spring(width=-20),
            Item("setting_panels", style="custom", show_label=False, width=-700),
            Spring(width=-20),
            orientation="horizontal"),
        kind="live",
        resizable=True)

    def __init__(self, figure, callback_dict):
        self.setting_panels = SettingPanels(callback_dict)
        self.figure_panel = FigurePanel(figure)

        self.edit_traits()
        app = QtGui.QApplication.instance()
        app = QtGui.QApplication.instance()
        for c in range(0,len(app.topLevelWidgets())):
            if hasattr(app.topLevelWidgets()[c], "_mw"):
                app.topLevelWidgets()[c]._mw.setMinimumSize(1200,600)
        start_event_loop_qt4()
