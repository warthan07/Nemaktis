#ifndef POSTPROCESSORSETTINGS_H
#define POSTPROCESSORSETTINGS_H

#include "BaseSettings.h"

class BasePostprocessorSettings : public BaseSettings {
public:
	BasePostprocessorSettings(
		const nlohmann::json &j, std::string section_name);

	const bool activate;
	const std::string base_name;
};

class MicrographOutputSettings : public BasePostprocessorSettings {
public:
	MicrographOutputSettings(const nlohmann::json &j);

	const std::vector<double> iso_layer_thickness;
	const std::vector<double> iso_layer_index;
	const double focalisation_z_shift;
	const double numerical_aperture;
};

class PostprocessorSettings : public BaseSettings {
public:
	PostprocessorSettings(const nlohmann::json &j);

	BasePostprocessorSettings volume_output;
	MicrographOutputSettings micrograph_output;
};

#endif
