#ifndef POSTPROCESSORSETTINGS_H
#define POSTPROCESSORSETTINGS_H

#include "BaseSettings.h"

class BasePostprocessorSettings : public BaseSettings {
public:
	BasePostprocessorSettings(
		const nlohmann::json &j, std::string section_name);

	/**
	 * Should we activate this postprocessor?
	 */
	const bool activate;
	/**
	 * Base name for the associated json section
	 */
	const std::string base_name;
};

class PostprocessorSettings : public BaseSettings {
public:
	PostprocessorSettings(const nlohmann::json &j);

	/**
	 * Json settings for volume output
	 */
	BasePostprocessorSettings volume_output;
	/**
	 * Json settings for screen output
	 */
	BasePostprocessorSettings micrograph_output;
};

#endif
