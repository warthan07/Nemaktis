#ifndef SCREENOUTPUT_H
#define SCREENOUTPUT_H

#include <Eigen/Core>

#include "PostprocessorBase.h"
#include "OpticalFieldCollection.h"

class ScreenOutput : public PostprocessorBase {
public:
	ScreenOutput(
		const RootSettings &settings,
		const PhysicsCoefficients &coefs);

	void apply(ScreenOpticalFieldCollection &screen_optical_fields);

	void apply_no_export(
		ScreenOpticalFieldCollection &screen_optical_fields,
		std::complex<double>* output_fields_vals);

private:
	/**
	 * Assemble the Fourier filter which propagate the optical fields
	 * forward through the isotropic layers and backward to the
	 * focalisation plane.
	 */
	std::shared_ptr<std::vector<Eigen::Matrix2cd> > assemble_fourier_filter(
		double wavelength, std::pair<double,double> q_val) const;

	const PhysicsCoefficients &coefs;

	int Nx, Ny, Nz;
	double delta_x, delta_y;

	std::vector<double> iso_layer_thickness, iso_layer_index;
	double z_foc, numerical_aperture;

	/**
	 * Array containing all the wavelengths in the light spectrum.
	 */
	std::vector<double> wavelengths;
	/**
	 * Array containing the incoming wavectors.
	 */
	std::vector<std::pair<double,double> > q_vals;
};
#endif
