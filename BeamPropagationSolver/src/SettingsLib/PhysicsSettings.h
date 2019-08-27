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

	const std::string ne_expression;
	const std::string no_expression;
	const std::string nhost_expression;
	const std::string nin_expression;
	std::vector<double> wavelengths() const {
		return _wavelengths;
	}

private:
	std::vector<double> _wavelengths;
};

class PhysicsSettings : public BaseSettings {
public:
	PhysicsSettings(const nlohmann::json &j);

	InitialConditionsSettings initial_conditions;
	CoefficientsSettings coefs;
};

#endif
