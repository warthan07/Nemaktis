#ifndef LCPROBLEMSETTINGS_H
#define LCPROBLEMSETTINGS_H

#include "BaseSettings.h"

enum class BasisConvention {
	XYZ, YZX, ZXY
};

class InitialConditionsSettings : public BaseSettings {
public:
	InitialConditionsSettings(const nlohmann::json &j);
	
	const std::string beam_profile_type;
	const std::string initial_solution_file;
	BasisConvention basis_convention() const {
		return _basis_convention;
	}

private:
	BasisConvention _basis_convention;
};

class CoefficientsSettings : public BaseSettings {
public:
	CoefficientsSettings(const nlohmann::json &j);

	const std::string no_expression;
	const std::string ne_expression;
	const std::string ne_imag_expression;
	const std::string nhost_expression;
	const std::string nin_expression;
	const std::string nout_expression;
	std::vector<double> wavelengths() const {
		return _wavelengths;
	}
	std::vector<std::pair<double,double> > q_vals() const {
		return _q_vals;
	}

private:
	std::vector<double> _wavelengths;
	std::vector<std::pair<double,double> > _q_vals;
};

class PhysicsSettings : public BaseSettings {
public:
	PhysicsSettings(const nlohmann::json &j);

	InitialConditionsSettings initial_conditions;
	CoefficientsSettings coefs;
};

#endif
