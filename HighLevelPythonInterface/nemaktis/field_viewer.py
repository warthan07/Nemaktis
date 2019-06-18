from .light_propagator import OpticalFields
from .field_viewer_ui import *
from dtmm import color as c

import os
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from scipy.interpolate import RectBivariateSpline

class FieldViewer:
    """A class allowing to recombine optical fields to generate optical micrographs like in a real
    microscope.

    Parameters
    ----------
    optical_fields : OpticalFields object
        Can be created either by a LightPropagator or directly by importing a vti file
        exported in a previous simulation.
    cmf : numpy ndarray
        A color matching array created with the `dtmm` package, 
        see `<https://dtmm.readthedocs.io/en/latest/reference.html#module-dtmm.color>`_

    """

    polariser = True 
    """Is there a polariser in the optical setup?"""
    analyser = True 
    """Is there an analyser in the optical setup?"""
    compensator = "No"
    """
    If "No", remove the compensator from the optical setup. Other values set the type of
    compensator:

    * "Quarter-wave": An achromatic quarter-wave compensator
    * "Half-wave": An achromatic half-wave compensator
    * "Tint-sensitive": a full-wave compensator at 540 nm.
    """

    polariser_angle = 0
    """Angle (in degree) between the privileged axis of the polariser and the x-axis"""
    analyser_angle = 90
    """Angle (in degree) between the privileged axis of the analyser and the x-axis"""
    compensator_angle = 0
    """Angle (in degree) between the fast axis of the compensator and the x-axis"""

    
    n_tiles_x = 1
    """Number of repetitions of the micrograph in the x-direction"""
    n_tiles_y = 1
    """Number of repetitions of the micrograph in the y-direction"""
    grayscale = False
    """Should we calculate a grayscale micrograph (True) or a color micrograph (False)"""

    intensity = 1
    """Intensity factor of the micrograph"""

    _slider_spans = {
        "intensity": [0, 2],
        "focus":     [-10, 10]}

    def __init__(self, optical_fields, cmf = None):
        if not isinstance(optical_fields, OpticalFields):
            raise TypeError("optical_fields should be an OpticalFields object")
        self._optical_fields = optical_fields.copy()

        shape = self._optical_fields.vals.shape
        self._image = np.zeros((shape[-2], shape[-1], 3))
        
        if cmf is None:
            self._cmf = c.load_tcmf(1000*optical_fields.get_wavelengths())
        else:
            self._cmf = cmf

    def plot(self):
        """Run the graphical user interface for real-time visualisation of micrographs"""

        print("{ Running field viewer graphical interface }")
        mpl.rcParams["toolbar"] = "None"
        self._fig = Figure(tight_layout=True)
        self._fig.set_facecolor("#efebe7")
        self._ax = self._fig.add_subplot(111)

        (self._Lx,self._Ly) = self._optical_fields.get_mesh_lengths()
        self._im_extent = (
            -self.n_tiles_x*self._Lx/2, self.n_tiles_x*self._Lx/2,
            -self.n_tiles_y*self._Ly/2, self.n_tiles_y*self._Ly/2)
    
        self.update_image()
        self._ax_im = self._ax.imshow(
            np.tile(self._image,(self.n_tiles_y,self.n_tiles_x,1)),
            origin="lower", interpolation="bicubic",
            extent=self._im_extent)

        def update_intensity(new_intensity):
            self.intensity = new_intensity
            self._update_plot()

        def update_focus(new_z):
            self._optical_fields.propagate(new_z)
            self._update_plot()

        def update_polariser(polariser):
            self.polariser = polariser=="Yes"
            self._update_plot()

        def update_polariser_angle(new_angle):
            self.polariser_angle = new_angle
            self._update_plot()

        def update_analyser(analyser):
            self.analyser = analyser=="Yes"
            self._update_plot()

        def update_analyser_angle(new_angle):
            self.analyser_angle = new_angle
            self._update_plot()

        def update_compensator(compensator):
            self.compensator = compensator
            self._update_plot()

        def update_compensator_angle(new_angle):
            self.compensator_angle = new_angle
            self._update_plot()

        def update_color_mode(color_mode):
            self.grayscale = color_mode=="Grayscale"
            self._update_plot()

        def update_n_tiles_x(new_int):
            self.n_tiles_x = new_int
            self._im_extent = (
                    -self.n_tiles_x*self._Lx/2, self.n_tiles_x*self._Lx/2,
                    -self.n_tiles_y*self._Ly/2, self.n_tiles_y*self._Ly/2)
            self._ax_im.set_extent(self._im_extent)
            self._update_plot()

        def update_n_tiles_y(new_int):
            self.n_tiles_y = new_int
            self._im_extent = (
                -self.n_tiles_x*self._Lx/2, self.n_tiles_x*self._Lx/2,
                -self.n_tiles_y*self._Ly/2, self.n_tiles_y*self._Ly/2)
            self._ax_im.set_extent(self._im_extent)
            self._update_plot()

        def save(file_path):
            im = np.tile(self._image,(self.n_tiles_y,self.n_tiles_x,1))
            interp_funcs = list(RectBivariateSpline(
                np.linspace(self._im_extent[2], self._im_extent[3], im.shape[0]),
                np.linspace(self._im_extent[0], self._im_extent[1], im.shape[1]),
                im[:,:,c], kx=3, ky=3) for c in range(0,3))
            dx = self._Lx*self.n_tiles_x/(im.shape[1]-1)
            dy = self._Ly*self.n_tiles_y/(im.shape[0]-1)
            new_Nx = int(self._Lx*self.n_tiles_x/min(dx,dy))+1
            new_Ny = int(self._Ly*self.n_tiles_y/min(dx,dy))+1
            new_ys = np.linspace(self._im_extent[2], self._im_extent[3], new_Ny)
            new_xs = np.linspace(self._im_extent[0], self._im_extent[1], new_Nx)
            new_im = np.concatenate(tuple(
                np.expand_dims(interp_funcs[c](new_ys,new_xs), 2) for c in range(0,3)), 2)
            plt.imsave(file_path, np.flip(np.clip(new_im, 0, 1), axis=0))

        ui = FieldViewerUI(self._fig, {
            "set_polariser": update_polariser,
            "set_polariser_angle": update_polariser_angle,
            "set_analyser": update_analyser,
            "set_analyser_angle": update_analyser_angle,
            "set_compensator": update_compensator,
            "set_compensator_angle": update_compensator_angle,
            "set_intensity": update_intensity,
            "set_z_focus": update_focus,
            "set_color_mode": update_color_mode,
            "set_n_tiles_x": update_n_tiles_x,
            "set_n_tiles_y": update_n_tiles_y,
            "save": save})



    @property
    def z_focus(self):
        """Current vertical position of the focal plane"""
        return self._optical_fields.z_focus


    @z_focus.setter
    def z_focus(self, new_z):
        self._optical_fields.propagate(new_z)


    def get_image(self):
        """Returns the current micrograph as a numpy array of shape (Ny,Nx,3|1),
        (last dim=3 if in color mode, 1 if in grayscale mode)."""
        return np.flip(np.tile(self._image,(self.n_tiles_y,self.n_tiles_x,1)), axis=0)


    def update_image(self):
        """Recompute the micrograph from the optical fields data"""
        fields_vals = self._optical_fields.vals

        if self.polariser:
            trans_mat = polariser_matrix(self.polariser_angle)
            trans_mat = np.matmul(sample_transfer_matrix(self._optical_fields.vals), trans_mat)
        else:
            trans_mat = sample_transfer_matrix(self._optical_fields.vals)

        if self.compensator=="Quarter-wave":
            trans_mat = np.matmul(quarter_waveplate_matrix(self.compensator_angle), trans_mat)
        elif self.compensator=="Half-wave":
            trans_mat = np.matmul(half_waveplate_matrix(self.compensator_angle), trans_mat)
        elif self.compensator=="Tint-sensitive":
            trans_mat = np.matmul(
                tint_sensitive_matrix(
                    self._optical_fields.get_wavelengths(),
                    self.compensator_angle),
                trans_mat)

        if self.analyser:
            trans_mat = np.matmul(polariser_matrix(self.analyser_angle), trans_mat)
        specter_data = 0.5*np.sum(
            np.abs(trans_mat)**2, axis=(-1,-2), keepdims=False)

        specter_data = self.intensity*specter_data.transpose((1,2,0))
        c.specter2color(specter_data, self._cmf, out = self._image, gray = self.grayscale)


    def _update_slider_pos(self):
        slider_pos = [[0.6, 0.71, 0.35, 0.03],
                      [0.6, 0.64, 0.35, 0.03],
                      [0.6, 0.57, 0.35, 0.03],
                      [0.6, 0.50, 0.35, 0.03]]
        if self.analyser:
            self._ax_widgets["analyser"].set_position(slider_pos.pop())
        if self.polariser:
            self._ax_widgets["polariser"].set_position(slider_pos.pop())
        self._ax_widgets["focus"].set_position(slider_pos.pop())
        self._ax_widgets["intensity"].set_position(slider_pos.pop())

        # to avoid overlap between widgets, we also reset the pos of the invisible sliders
        if not self.analyser:
            self._ax_widgets["analyser"].set_position(slider_pos.pop())
        if not self.polariser:
            self._ax_widgets["polariser"].set_position(slider_pos.pop())


    def _update_slider_span(self, name):
        valmin = self._slider_spans[name][0]
        valmax = self._slider_spans[name][1]
        self._widgets[name].valmin = valmin
        self._widgets[name].valmax = valmax
        self._ax_widgets[name].set_xlim(valmin, valmax)
        self._fig.canvas.draw_idle()


    def _update_plot(self):
        self.update_image()
        self._ax_im.set_data(np.tile(self._image,(self.n_tiles_y,self.n_tiles_x,1)))
        self._fig.canvas.draw_idle()


def polariser_matrix(angle):
    cosA = np.cos(angle*np.pi/180)
    sinA = np.sin(angle*np.pi/180)
    return np.array([
        [cosA**2, cosA*sinA],
        [cosA*sinA, sinA**2]])[np.newaxis, np.newaxis, np.newaxis, :, :]

def quarter_waveplate_matrix(angle):
    cosA = np.cos(angle*np.pi/180)
    sinA = np.sin(angle*np.pi/180)
    return np.exp(-1j*np.pi/4)*np.array([
        [cosA**2+1j*sinA**2, (1-1j)*cosA*sinA],
        [(1-1j)*cosA*sinA, sinA**2+1j*cosA**2]])[np.newaxis, np.newaxis, np.newaxis, :, :]

def half_waveplate_matrix(angle):
    cosA = np.cos(angle*np.pi/180)
    sinA = np.sin(angle*np.pi/180)
    return np.exp(-1j*np.pi/2)*np.array([
        [cosA**2-sinA**2, 2*cosA*sinA],
        [2*cosA*sinA, sinA**2-cosA**2]])[np.newaxis, np.newaxis, np.newaxis, :, :]

def tint_sensitive_matrix(wavelengths, angle):
    cosA = np.cos(angle*np.pi/180)
    sinA = np.sin(angle*np.pi/180)
    expG = np.exp(1j*np.pi*0.540/wavelengths) # Full-wave at 540 nm
    return np.array([
        [np.conj(expG)*cosA**2+expG*sinA**2, (np.conj(expG)-expG)*cosA*sinA],
        [(np.conj(expG)-expG)*cosA*sinA, np.conj(expG)*sinA**2+expG*cosA**2]]).transpose(
            (2,0,1))[:, np.newaxis, np.newaxis, :, :]

def sample_transfer_matrix(optical_vals):
    dims = optical_vals.shape
    return optical_vals.reshape((dims[0],2,2,dims[2],dims[3])).transpose((0,3,4,2,1))
