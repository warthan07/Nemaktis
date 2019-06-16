#ifndef BEAMPROFILE_H
#define BEAMPROFILE_H

#include <complex>

#include "PhysicsCoefficients.h"

enum class BeamProfileType {
	GaussianBeam,
	UniformBeam
};

class BeamProfile {
public:
	BeamProfile(const PhysicsCoefficients &coefs, double pol_angle, double wavelength);

	virtual std::complex<double> get_Ex(double x, double y) = 0;
	virtual std::complex<double> get_Ey(double x, double y) = 0;

protected:
	/**
	 * Extraordinary refractive index
	 */
	double ne;
	/**
	 * Ordinary refractive index
	 */
	double no;
	/**
	 * Wavevector in empty space
	 */
	double k;
	/**
	 * Input polarisation 
	 */
	double theta;
	/**
	 * Sample thickness
	 */
	double h;
};

class GaussianBeam : public BeamProfile {
public:
	GaussianBeam(
		const PhysicsCoefficients &coefs, double waist, 
		double pol_angle, double wavelength);

	virtual std::complex<double> get_Ex(double x, double y);
	virtual std::complex<double> get_Ey(double x, double y);

private:
	/**
	 * Beam waist
	 */
	const double w;
};

class UniformBeam : public BeamProfile {
public:
	UniformBeam(const PhysicsCoefficients &coefs, double pol_angle, double wavelength);

	virtual std::complex<double> get_Ex(double x, double y);
	virtual std::complex<double> get_Ey(double x, double y);
};

#endif
