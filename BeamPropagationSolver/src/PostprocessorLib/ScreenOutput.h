#ifndef SCREENOUTPUT_H
#define SCREENOUTPUT_H

#include "PostprocessorBase.h"

class ScreenOutput : public PostprocessorBase {
public:
	ScreenOutput(
		const RootSettings &settings,
		const PhysicsCoefficients &coefs);

	virtual void apply(
		VectorField<double> &lc_sol,
		std::vector<VectorField<std::complex<double> > > (&bpm_sol)[2]);

	void apply_no_export(
		VectorField<double> &lc_sol,
		std::vector<VectorField<std::complex<double> > > (&bpm_sol)[2],
		std::complex<double>* fields_vals);

private:
	/**
	 * Assemble the Fourier filter which propagate the optical fields
	 * forward through the isotropic layers and backward to the
	 * focalisation plane.
	 */
	std::shared_ptr<std::vector<std::complex<double> > > assemble_iso_filter(
		double wavelength) const;

	int Nx, Ny, Nz;
	double delta_x, delta_y;

	std::vector<double> iso_layer_thickness;
	std::vector<double> iso_layer_index;
	double focalisation_z_shift;
	double numerical_aperture;

	std::vector<double> wavelengths;
};
#endif
