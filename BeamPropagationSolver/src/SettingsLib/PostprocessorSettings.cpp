#include "PostprocessorSettings.h"

BasePostprocessorSettings::BasePostprocessorSettings(
		const nlohmann::json &j, std::string section_name) :
	BaseSettings(j, section_name),
	activate(parse<bool>("Activate")),
	base_name(parse<std::string>("Base name")) {}

PostprocessorSettings::PostprocessorSettings(const nlohmann::json &j) :
	BaseSettings(j, "Postprocessor settings"),
	volume_output(json_node, "Bulk output"),
	micrograph_output(json_node, "Screen output") {}
