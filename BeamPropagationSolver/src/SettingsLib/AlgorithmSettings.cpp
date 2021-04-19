#include "AlgorithmSettings.h"

GeneralSettings::GeneralSettings(const nlohmann::json &j) :
	BaseSettings(j, "General"),
	optical_axis_field_type(parse<std::string>("Optical axis field type")),
	results_folder_name(parse<std::string>("Results folder name")) {

	if(optical_axis_field_type!="Director" && optical_axis_field_type!="Q-tensor" &&
			optical_axis_field_type!="DualDirector" && optical_axis_field_type!="None")
		throw std::string(
			"Parsing error: \"LC field type\" must be \"Director\", \"Q-tensor\", "
			"\"DualDirector\" or \"None\"");
}

BPMSettings::BPMSettings(const nlohmann::json &j) :
	BaseSettings(j, "Beam propagation"),
	N_woodbury_steps(parse<int>("N Woodbury steps")),
	Nz_substeps(parse<int>("Number of substeps per slab")){}


AlgorithmSettings::AlgorithmSettings(const nlohmann::json &j) :
	BaseSettings(j, "Algorithm settings"),
	general(json_node),
	bpm(json_node) {}
