/** 
 * @file PhysicsCoefficients.h
 */

#ifndef LCCOEFFICIENTS_H
#define LCCOEFFICIENTS_H

#include <memory>
#include <map>

#include "RootSettings.h"
#include "CartesianMesh.h"

#define PI 3.1415926535897932

class PhysicsCoefficients {
public:
	PhysicsCoefficients(
		const RootSettings &settings,
		const CartesianMesh &mesh);

	/**
	 * Principal refractive index of selected sample domain at the specified wavelength.
	 *
	 * Relevant for isotropic, uniaxial and biaxial domains.
	 */
	double get_n1(unsigned int domain_id, double wavelength) const;
	/**
	 * Secondary refractive index of selected sample domain at the specified wavelength.
	 *
	 * Relevant for uniaxial and biaxial domains.
	 */
	double get_n2(unsigned int domain_id, double wavelength) const;
	/**
	 * Tertiary refractive index of selected sample domain at the specified wavelength.
	 *
	 * Relevant for biaxial domains.
	 */
	double get_n3(unsigned int domain_id, double wavelength) const;

	/**
	 * Refractive index of selected upper isotropic layer at the specified wavelength.
	 */
	double get_niso_up(unsigned int layer_id, double wavelength) const;
	/**
	 * Refractive index of selected lower isotropic layer at the specified wavelength.
	 */
	double get_niso_lo(unsigned int layer_id, double wavelength) const;

	/**
	 * Refractive index of input medium at the specified wavelength
	 */
	double get_nin(double wavelength) const;
	/**
	 * Refractive index of input medium at the specified wavelength
	 */
	double get_nout(double wavelength) const;

	/**
	 * Wavelengths of input plane waves
	 */
	const std::vector<double> &wavelengths() const {
		return coefs_settings.wavelengths();
	}
	/**
	 * Transverse wavevevectors of input plane waves renormalized by k0=2pi/wavelength
	 */
	const std::vector<std::pair<double,double> > &q_vals() const {
		return coefs_settings.q_vals();
	}
	/**
	 * Thicknesses of the upper isotropic layer
	 */
	const std::vector<double> &hiso_up() const {
		return coefs_settings.hiso_up();
	}
	/**
	 * Thicknesses of the lower isotropic layer
	 */
	const std::vector<double> &hiso_lo() const {
		return coefs_settings.hiso_lo();
	}

private:
	const CoefficientsSettings &coefs_settings;
};

#endif
