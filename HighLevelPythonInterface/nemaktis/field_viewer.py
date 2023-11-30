from .light_propagator import OpticalFields
from .field_viewer_ui import *
from dtmm import color

import os
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from scipy.interpolate import RectBivariateSpline

import time

class FieldViewer:
    """A class allowing to recombine optical fields to generate optical micrographs like in
    a real microscope. For more details, see `[Imaging of the object] 
    <https://nemaktis.readthedocs.io/en/latest/intro/microscopy_model.html#imaging-of-the-object>`_
    and `[Optical elements for polarized optical micrographs]
    <https://nemaktis.readthedocs.io/en/latest/intro/microscopy_model.html#optical-elements-for-polarized-optical-micrographs>`_

    Parameters
    ----------
    optical_fields : :class:`~nemaktis.light_propagator.OpticalFields` object
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
    upper_waveplate = "No"
    """
    If "No", remove the upper waveplate from the optical setup. Other values set the type of
    waveplate:

    * "Quarter-wave": An achromatic quarter-wave compensator
    * "Half-wave": An achromatic half-wave compensator
    * "Tint-sensitive": a full-wave compensator at 540 nm.
    """
    lower_waveplate = "No"
    """
    If "No", remove the upper waveplate from the optical setup. Other values set the type of
    waveplate:

    * "Quarter-wave": An achromatic quarter-wave compensator
    * "Half-wave": An achromatic half-wave compensator
    * "Tint-sensitive": a full-wave compensator at 540 nm.
    """

    polariser_angle = 0
    """Angle (in degree) between the privileged axis of the polariser and the x-axis"""
    analyser_angle = 90
    """Angle (in degree) between the privileged axis of the analyser and the x-axis"""
    upper_waveplate_angle = 0
    """Angle (in degree) between the fast axis of the upper waveplate and the x-axis"""
    lower_waveplate_angle = 0
    """Angle (in degree) between the fast axis of the lower waveplate and the x-axis"""
    angle_lock = False
    """Should the relative angles between optical elements be locked?"""

    intensity = 1
    """Intensity factor of the micrograph"""
    NA_condenser = 0
    """Numerical aperture of the microscope's condenser"""
    
    n_tiles_x = 1
    """Number of repetitions of the micrograph in the x-direction"""
    n_tiles_y = 1
    """Number of repetitions of the micrograph in the y-direction"""
    grayscale = False
    """Should we calculate a grayscale micrograph (True) or a color micrograph (False)"""


    def __init__(self, optical_fields, cmf = None):
        if not isinstance(optical_fields, OpticalFields):
            raise TypeError("optical_fields should be an OpticalFields object")
        self._optical_fields = optical_fields
        optical_fields.focus_fields(0)
        self._NA_objective = optical_fields._max_NA_objective

        shape = self._optical_fields.vals.shape
        self._specter = np.zeros((shape[3],shape[4],shape[0],shape[1]))
        self._q_averaged_specter = np.zeros((shape[3],shape[4],shape[0]))
        self._image = np.zeros((shape[3], shape[4], 3))
        
        if cmf is None:
            self._cmf = color.load_tcmf(1000*optical_fields.get_wavelengths())
        else:
            self._cmf = cmf

        self._koehler_1D = self._optical_fields._koehler_1D
        self._N_radial_wavevectors = self._optical_fields._N_radial_wavevectors

    def plot(self):
        """Run a graphical user interface allowing to dynamically adjust
        the attributes of this class and visualize the associated
        micrographs in real-time."""

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

        def update_all_angles(angle_increment):
            new_angles = {}
            if self.polariser:
                self.polariser_angle = (self.polariser_angle+angle_increment+90)%180-90
                new_angles["polariser"] = self.polariser_angle
            if self.analyser:
                self.analyser_angle = (self.analyser_angle+angle_increment+90)%180-90
                new_angles["analyser"] = self.analyser_angle
            if self.lower_waveplate!="No":
                self.lower_waveplate_angle = (self.lower_waveplate_angle+angle_increment+90)%180-90
                new_angles["lower_waveplate"] = self.lower_waveplate_angle
            if self.upper_waveplate!="No":
                self.upper_waveplate_angle = (self.upper_waveplate_angle+angle_increment+90)%180-90
                new_angles["upper_waveplate"] = self.upper_waveplate_angle
            return new_angles

        def update_polariser(polariser):
            self.polariser = polariser=="Yes"
            self.update_image()
            self._update_plot()

        def update_polariser_angle(new_angle):
            if self.angle_lock:
                new_angles = update_all_angles(new_angle-self.polariser_angle)
            else:
                self.polariser_angle = new_angle
                new_angles = {}
            self.update_image()
            self._update_plot()
            return new_angles

        def update_analyser(analyser):
            self.analyser = analyser=="Yes"
            self.update_image()
            self._update_plot()

        def update_analyser_angle(new_angle):
            if self.angle_lock:
                new_angles = update_all_angles(new_angle-self.analyser_angle)
            else:
                self.analyser_angle = new_angle
                new_angles = {}
            self.update_image()
            self._update_plot()
            return new_angles

        def update_upper_waveplate(waveplate):
            self.upper_waveplate = waveplate
            self.update_image()
            self._update_plot()

        def update_upper_waveplate_angle(new_angle):
            if self.angle_lock:
                new_angles = update_all_angles(new_angle-self.upper_waveplate_angle)
            else:
                self.upper_waveplate_angle = new_angle
                new_angles = {}
            self.update_image()
            self._update_plot()
            return new_angles

        def update_lower_waveplate(waveplate):
            self.lower_waveplate = waveplate
            self.update_image()
            self._update_plot()

        def update_lower_waveplate_angle(new_angle):
            if self.angle_lock:
                new_angles = update_all_angles(new_angle-self.lower_waveplate_angle)
            else:
                self.lower_waveplate_angle = new_angle
                new_angles = {}
            self.update_image()
            self._update_plot()
            return new_angles

        def update_angle_lock(angle_lock):
            self.angle_lock = angle_lock=="Yes"

        def update_intensity(new_intensity):
            self.intensity = new_intensity
            self._update_POM_image()
            self._update_plot()

        def update_focus(new_z):
            self.z_focus = new_z
            self.update_image()
            self._update_plot()

        def update_NA_condenser(new_NA):
            if new_NA>self.NA_condenser:
                self.NA_condenser = max(0,min(self._optical_fields._max_NA_condenser,new_NA))
                self._update_specter()
            else:
                self.NA_condenser = max(0,min(self._optical_fields._max_NA_condenser,new_NA))
            self._update_averaged_specter()
            self._update_POM_image()
            self._update_plot()

        def update_NA_objective(new_NA):
            self._NA_objective = new_NA
            self._optical_fields.update_NA_objective(new_NA)
            self._optical_fields.focus_fields()
            self.update_image()
            self._update_plot()

        def update_color_mode(color_mode):
            self.grayscale = color_mode=="Grayscale"
            self._update_POM_image()
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

        ui = FieldViewerUI(
            self._fig, {
            "set_polariser": update_polariser,
            "set_polariser_angle": update_polariser_angle,
            "set_lower_waveplate": update_lower_waveplate,
            "set_lower_waveplate_angle": update_lower_waveplate_angle,
            "set_upper_waveplate": update_upper_waveplate,
            "set_upper_waveplate_angle": update_upper_waveplate_angle,
            "set_analyser": update_analyser,
            "set_analyser_angle": update_analyser_angle,
            "set_angle_lock": update_angle_lock,
            "set_intensity": update_intensity,
            "set_z_focus": update_focus,
            "set_NA_condenser": update_NA_condenser,
            "set_NA_objective": update_NA_objective,
            "set_color_mode": update_color_mode,
            "set_n_tiles_x": update_n_tiles_x,
            "set_n_tiles_y": update_n_tiles_y,
            "save": save}, {
            "polariser": self.polariser,
            "lower_waveplate": self.lower_waveplate,
            "upper_waveplate": self.upper_waveplate,
            "analyser": self.analyser,
            "polariser_angle": self.polariser_angle,
            "lower_waveplate_angle": self.lower_waveplate_angle,
            "upper_waveplate_angle": self.upper_waveplate_angle,
            "analyser_angle": self.analyser_angle,
            "angle_lock": self.angle_lock,
            "intensity": self.intensity,
            "z_focus": self.z_focus,
            "max_NA_condenser": self._optical_fields._max_NA_condenser,
            "max_NA_objective": self._optical_fields._max_NA_objective,
            "NA_condenser": self.NA_condenser,
            "NA_objective": self._NA_objective,
            "n_tiles_x": self.n_tiles_x,
            "n_tiles_y": self.n_tiles_y,
            "grayscale": self.grayscale})


    @property
    def z_focus(self):
        """Current vertical position of the focal plane"""
        return self._optical_fields.z_focus


    @z_focus.setter
    def z_focus(self, new_z):
        self._optical_fields.focus_fields(new_z)


    @property
    def NA_objective(self):
        """Numerical aperture of the microscope's objective"""
        return self._NA_objective


    @z_focus.setter
    def NA_objective(self, new_NA):
        self._NA_objective = new_NA
        self._optical_fields.update_NA_objective(new_NA)
        self._optical_fields.focus_fields()


    def get_image(self):
        """Returns the current micrograph as a numpy array of shape (Ny,Nx,3|1),
        (last dim is 3 if in color mode, 1 if in grayscale mode)."""
        return np.flip(np.tile(self._image,(self.n_tiles_y,self.n_tiles_x,1)), axis=0)

    def get_spectrum(self):
        """Returns the q-averaged spectrum as a numpy array of shape (Ny,Nx,Nl),
        (last dim is 3 if in color mode, 1 if in grayscale mode)."""
        return np.flip(np.tile(self._q_averaged_specter,(self.n_tiles_y,self.n_tiles_x,1)), axis=0)


    def update_image(self):
        """Recompute the micrograph from the optical fields data"""
        self._update_specter()
        self._update_averaged_specter()
        self._update_POM_image()


    def _update_specter(self):
        qr_idx_max = min(
            self._N_radial_wavevectors-1,
            self._optical_fields.get_qr_index(self.NA_condenser))
        if not self._koehler_1D:
            q_idx_start = 0
            q_idx_end = get_q_idx(qr_idx_max+1)
        else:
            q_idx_start = self._N_radial_wavevectors-1-qr_idx_max
            q_idx_end = self._N_radial_wavevectors+qr_idx_max
        fields_vals = self._optical_fields.focused_vals[:,q_idx_start:q_idx_end,:,:,:]

        if not self.polariser and self.lower_waveplate=="No" and self.upper_waveplate=="No" and not self.analyser:
            self._specter[:,:,:,q_idx_start:q_idx_end] = \
                0.5*np.sum(np.abs(fields_vals)**2, axis=2, keepdims=False).transpose((2,3,0,1))
        else:
            if self.polariser:
                trans_mat = polariser_matrix(self.polariser_angle)
                if self.lower_waveplate!="No":
                    trans_mat = np.matmul(waveplate_matrix(
                        self.lower_waveplate, self.lower_waveplate_angle,
                        self._optical_fields.get_wavelengths()), trans_mat)
                trans_mat = np.matmul(sample_transfer_matrix(fields_vals), trans_mat)
            elif self.lower_waveplate!="No":
                trans_mat = waveplate_matrix(
                    self.lower_waveplate, self.lower_waveplate_angle,
                    self._optical_fields.get_wavelengths())
                trans_mat = np.matmul(sample_transfer_matrix(fields_vals), trans_mat)
            else:
                trans_mat = sample_transfer_matrix(fields_vals)

            if self.upper_waveplate!="No":
                trans_mat = np.matmul(waveplate_matrix(
                    self.upper_waveplate, self.upper_waveplate_angle,
                    self._optical_fields.get_wavelengths()), trans_mat)

            if self.analyser:
                trans_mat = np.matmul(polariser_matrix(self.analyser_angle), trans_mat)
            specter_data = 0.5*np.sum(np.abs(trans_mat)**2, axis=(-1,-2), keepdims=False)

            self._specter[:,:,:,q_idx_start:q_idx_end] = specter_data.transpose((2,3,0,1))


    def _update_averaged_specter(self):
        if self.NA_condenser<1e-8:
            # For vanishingly small NA, we only need the data from the wavector q=(0,0)
            # It is necessary to distinguish this case since we renormalize the integration
            # weights by the NA of the condenser in the following general case
            q_idx = self._N_radial_wavevectors-1 if self._koehler_1D else 0
            self._q_averaged_specter = self._specter[:,:,:,q_idx]
        else:
            # We assemble the integration weights for the integral:
            #   \int specter(q,theta)*g(q) dq d(theta), with g(q)=q/(1-q**2)
            # specter is linearly interpolated along q and averaged between [0,2pi] along
            # theta. The maximum value for q is simply the NA of the condenser.  We use
            # exact primitives related to g(q) to get an acurate integration even when the
            # NA of the condenser is smoothly varied
            qr_idx_max = min(
                self._N_radial_wavevectors-1,
                self._optical_fields.get_qr_index(self.NA_condenser))
            delta_qr = self._optical_fields.get_delta_qr()
            if not self._koehler_1D:
                integration_weights = np.zeros((get_q_idx(qr_idx_max+1)))
                if qr_idx_max==1:
                    integration_weights[0] = \
                        (q_int_func(self.NA_condenser,delta_qr)-q_int_func(0,delta_qr))/delta_qr
                    integration_weights[1:7] = \
                        (q_int_func(0,0)-q_int_func(self.NA_condenser,0))/(6*delta_qr)
                else:
                    integration_weights[0] = \
                        (q_int_func(delta_qr,delta_qr)-q_int_func(0,delta_qr))/delta_qr
                    for qr_idx in range(1,qr_idx_max-1):
                        integration_weights[get_q_idx(qr_idx):get_q_idx(qr_idx+1)] = \
                            ( q_int_func((qr_idx+1)*delta_qr,(qr_idx+1)*delta_qr) \
                            + q_int_func((qr_idx-1)*delta_qr,(qr_idx-1)*delta_qr) \
                            - 2*q_int_func(qr_idx*delta_qr,qr_idx*delta_qr) ) \
                                / (6*qr_idx*delta_qr)
                    integration_weights[get_q_idx(qr_idx_max-1):get_q_idx(qr_idx_max)] = \
                        ( q_int_func(self.NA_condenser,qr_idx_max*delta_qr) \
                        + q_int_func((qr_idx_max-2)*delta_qr,(qr_idx_max-2)*delta_qr) \
                        - 2*q_int_func((qr_idx_max-1)*delta_qr,(qr_idx_max-1)*delta_qr) ) \
                            / (6*(qr_idx_max-1)*delta_qr)
                    integration_weights[get_q_idx(qr_idx_max):get_q_idx(qr_idx_max+1)] = \
                        ( q_int_func((qr_idx_max-1)*delta_qr,(qr_idx_max-1)*delta_qr) \
                        - q_int_func(self.NA_condenser,(qr_idx_max-1)*delta_qr) ) \
                            / (6*qr_idx_max*delta_qr)
                self._q_averaged_specter = np.sum(
                    self._specter[:,:,:,0:get_q_idx(qr_idx_max+1)] * (2/self.NA_condenser**2) *
                    integration_weights[np.newaxis,np.newaxis,np.newaxis,:],
                    axis=3, keepdims=False)
            else:
                Nr = self._N_radial_wavevectors
                delta_qr_edge = self.NA_condenser-(qr_idx_max-1)*delta_qr

                integration_weights = np.zeros(2*qr_idx_max+1)
                integration_weights[:] = delta_qr
                if qr_idx_max==1:
                    integration_weights[[0,-1]] = delta_qr_edge**2/(2*delta_qr)
                    integration_weights[1] = delta_qr_edge*( 2-delta_qr_edge/delta_qr )
                else:
                    integration_weights[[0,-1]] = delta_qr_edge**2/(2*delta_qr)
                    integration_weights[[1,-2]] = delta_qr_edge*( 1-delta_qr_edge/(2*delta_qr) ) + delta_qr/2

                q_vals = np.arange(-qr_idx_max,qr_idx_max+1)*delta_qr
                self._q_averaged_specter = np.sum(
                    self._specter[:,:,:,(Nr-1-qr_idx_max):(Nr+qr_idx_max)] / (2*self.NA_condenser) *
                    integration_weights[np.newaxis,np.newaxis,np.newaxis,:],
                    axis=3, keepdims=False)


    def _update_POM_image(self):
        color.specter2color(
            self.intensity*self._q_averaged_specter, self._cmf,
            out = self._image, gray = self.grayscale)
 

    def _update_plot(self):
        if self._fig.canvas is not None:
            self._ax_im.set_data(np.tile(self._image,(self.n_tiles_y,self.n_tiles_x,1)))
            self._fig.canvas.draw_idle()


def polariser_matrix(angle):
    cosA = np.cos(angle*np.pi/180)
    sinA = np.sin(angle*np.pi/180)
    return np.array([
        [cosA**2, cosA*sinA],
        [cosA*sinA, sinA**2]])[np.newaxis,np.newaxis,np.newaxis,np.newaxis,:,:]

def quarter_waveplate_matrix(angle):
    cosA = np.cos(angle*np.pi/180)
    sinA = np.sin(angle*np.pi/180)
    return np.exp(-1j*np.pi/4)*np.array([
        [cosA**2+1j*sinA**2, (1-1j)*cosA*sinA],
        [(1-1j)*cosA*sinA, sinA**2+1j*cosA**2]])[np.newaxis,np.newaxis,np.newaxis,np.newaxis,:,:]

def half_waveplate_matrix(angle):
    cosA = np.cos(angle*np.pi/180)
    sinA = np.sin(angle*np.pi/180)
    return np.exp(-1j*np.pi/2)*np.array([
        [cosA**2-sinA**2, 2*cosA*sinA],
        [2*cosA*sinA, sinA**2-cosA**2]])[np.newaxis,np.newaxis,np.newaxis,np.newaxis,:,:]

def tint_sensitive_matrix(wavelengths, angle):
    cosA = np.cos(angle*np.pi/180)
    sinA = np.sin(angle*np.pi/180)
    expG = np.exp(1j*np.pi*0.540/wavelengths) # Full-wave at 540 nm
    #  expG = np.exp(1j*np.pi*(1.25*(0.540/wavelengths-1)+1)) # Andrew's full-wave plate
    return np.array([
        [np.conj(expG)*cosA**2+expG*sinA**2, (np.conj(expG)-expG)*cosA*sinA],
        [(np.conj(expG)-expG)*cosA*sinA, np.conj(expG)*sinA**2+expG*cosA**2]]).transpose(
            (2,0,1))[:,np.newaxis,np.newaxis,np.newaxis,:,:]

def waveplate_matrix(waveplate_type, angle, wavelengths):
    if waveplate_type=="Quarter-wave":
        return quarter_waveplate_matrix(angle)
    elif waveplate_type=="Half-wave":
        return half_waveplate_matrix(angle)
    elif waveplate_type=="Tint-sensitive":
        return tint_sensitive_matrix(wavelengths, angle)

def sample_transfer_matrix(optical_vals):
    dims = optical_vals.shape
    return optical_vals.reshape((dims[0],dims[1],2,2,dims[3],dims[4])).transpose((0,1,4,5,3,2))

def q_int_func(q,qp):
    return q-np.arctanh(q)-qp*np.log(1-q*q)/2
def get_q_idx(qr_idx):
    return 1+3*qr_idx*(qr_idx-1)
