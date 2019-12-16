#include "AlgorithmSettings.h"

GeneralSettings::GeneralSettings(const nlohmann::json &j) :
	BaseSettings(j, "General"),
	lc_field_type(parse<std::string>("LC field type")),
	results_folder_name(parse<std::string>("Results folder name")) {

	if(lc_field_type!="Director" && lc_field_type!="Q-tensor")
		throw std::string(
			"Parsing error: \"LC field type\" must be \"Director\" or \"Q-tensor\"");
}

BPMSettings::BPMSettings(const nlohmann::json &j) :
	BaseSettings(j, "Beam propagation"),
	n_woodbury_steps(parse<int>("N Woodbury steps")),
	Nz_substeps(parse<int>("Number of substeps per slab")){}


AlgorithmSettings::AlgorithmSettings(const nlohmann::json &j) :
	BaseSettings(j, "Algorithm settings"),
	general(json_node),
	bpm(json_node) {}
