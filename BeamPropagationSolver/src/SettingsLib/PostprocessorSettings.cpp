#include "PostprocessorSettings.h"

BasePostprocessorSettings::BasePostprocessorSettings(
		const nlohmann::json &j, std::string section_name) :
	BaseSettings(j, section_name),
	activate(parse<bool>("Activate")),
	base_name(parse<std::string>("Base name")) {}

MicrographOutputSettings::MicrographOutputSettings(const nlohmann::json &j) :
	BasePostprocessorSettings(j, "Screen output"),
	iso_layer_thickness(parse_vector<double>("Isotropic layer thicknesses", 0)),
	iso_layer_index(parse_vector<double>("Isotropic layer refractive indices", 0)),
	focalisation_z_shift(parse<double>("Focalisation z-shift")),
	numerical_aperture(parse<double>("Numerical aperture")) {}

PostprocessorSettings::PostprocessorSettings(const nlohmann::json &j) :
	BaseSettings(j, "Postprocessor settings"),
	volume_output(json_node, "Bulk output"),
	micrograph_output(json_node) {}
