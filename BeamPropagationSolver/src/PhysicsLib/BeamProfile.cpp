#include "BeamProfile.h"

BeamProfile::BeamProfile(
		const PhysicsCoefficients &coefs, double pol_angle, double wavelength) :
	ne(coefs.get_ne(wavelength)),
	no(coefs.get_no(wavelength)),
	k(2*PI/wavelength),
	theta(pol_angle),
	h(coefs.mesh_thickness()) {}

GaussianBeam::GaussianBeam(
		const PhysicsCoefficients &coefs, double waist,
		double pol_angle, double wavelength) :
	BeamProfile(coefs, pol_angle, wavelength),
	w(waist) {}

std::complex<double> GaussianBeam::get_Ex(double x, double y) {

	return std::cos(theta) * std::exp(-(x*x+y*y)/(2*w*w));
}

std::complex<double> GaussianBeam::get_Ey(double x, double y) {

	return std::sin(theta) * std::exp(-(x*x+y*y)/(2*w*w));
}

UniformBeam::UniformBeam(
		const PhysicsCoefficients &coefs, double pol_angle, double wavelength) :
	BeamProfile(coefs, pol_angle, wavelength) {}

std::complex<double> UniformBeam::get_Ex(double x, double y) {

	return std::cos(theta)*std::exp(std::complex<double>(0,1)*0.*k*y);
}

std::complex<double> UniformBeam::get_Ey(double x, double y) {

	return std::sin(theta)*std::exp(std::complex<double>(0,1)*0.*k*y);
}
