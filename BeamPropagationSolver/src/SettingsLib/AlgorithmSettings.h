#ifndef ALGORITHMSETTINGS_H
#define ALGORITHMSETTINGS_H

#include "BaseSettings.h"

enum class BoundaryType {
	Transparent, Periodic
};

class GeneralSettings : public BaseSettings {
public:
	GeneralSettings(const nlohmann::json &j);

	const std::string lc_field_type;
	const std::string results_folder_name;
};

class BPMSettings : public BaseSettings {
public:
	BPMSettings(const nlohmann::json &j);

	std::vector<bool> wide_angle_corrections() const {
		return _wide_angle_corrections;
	}
	std::vector<BoundaryType> boundary_types() const {
		return _boundary_types;
	}

	const unsigned int n_woodbury_steps;

private:
	std::vector<bool> _wide_angle_corrections;
	std::vector<BoundaryType> _boundary_types;
};

class AlgorithmSettings : public BaseSettings {
public:
	AlgorithmSettings(const nlohmann::json &j);

	GeneralSettings general;
	BPMSettings bpm;
};

#endif
