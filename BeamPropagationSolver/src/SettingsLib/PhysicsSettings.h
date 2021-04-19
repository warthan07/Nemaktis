#ifndef LCPROBLEMSETTINGS_H
#define LCPROBLEMSETTINGS_H

#include "BaseSettings.h"

enum class BasisConvention {
	XYZ, YZX, ZXY
};

class FieldSetupSettings : public BaseSettings {
public:
	FieldSetupSettings(const nlohmann::json &j);
	
	/**
	 * Type of input beam profile
	 */
	const std::string beam_profile_type;
	/**
	 * Path to the microscope sample input vti file
	 */
	const std::string sample_file;

	/**
	 * Orientation convention for the microscope sample file
	 */
	BasisConvention basis_convention() const {
		return _basis_convention;
	}

private:
	BasisConvention _basis_convention;
};

class CoefficientsSettings : public BaseSettings {
public:
	CoefficientsSettings(
		const nlohmann::json &j,
		const std::string &sample_file);

	/**
	 * Wavelengths of incoming plane waves
	 */
	const std::vector<double> &wavelengths() const {
		return _wavelengths;
	}
	/**
	 * Transverse wavevectors of incoming plane waves renormalized by k0=2pi/wavelength
	 */
	const std::vector<std::pair<double,double> > &q_vals() const {
		return _q_vals;
	}
	/**
	 * Expressions of the principal refractive indices of all domains inside the sample
	 */
	const std::vector<std::string> &nsample_expressions() const {
		return _nsample_expressions;
	}
	/**
	 * Expressions of the refractive indices of the upper isotropic layer
	 */
	const std::vector<std::string> &niso_up_expressions() const {
		return _niso_up_expressions;
	}
	/**
	 * Expressions of the refractive indices of the lower isotropic layer
	 */
	const std::vector<std::string> &niso_lo_expressions() const {
		return _niso_lo_expressions;
	}
	/**
	 * Thicknesses of the upper isotropic layer
	 */
	const std::vector<double> &hiso_up() const {
		return _hiso_up;
	}
	/**
	 * Thicknesses of the lower isotropic layer
	 */
	const std::vector<double> &hiso_lo() const {
		return _hiso_lo;
	}
	/**
	 * Expression of the refractive index of the input medium
	 */
	const std::string nin_expression;
	/**
	 * Expression of the refractive index of the output medium
	 */
	const std::string nout_expression;

private:
	std::vector<std::string> _nsample_expressions;
	std::vector<std::string> _niso_up_expressions;
	std::vector<std::string> _niso_lo_expressions;
	std::vector<double> _hiso_up;
	std::vector<double> _hiso_lo;
	std::vector<double> _wavelengths;
	std::vector<std::pair<double,double> > _q_vals;
};

class PhysicsSettings : public BaseSettings {
public:
	PhysicsSettings(const nlohmann::json &j);

	FieldSetupSettings field_setup;
	CoefficientsSettings coefs;
};

#endif
