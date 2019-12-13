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
	n_woodbury_steps(parse<unsigned int>("N Woodbury steps")) {

	auto types = parse_vector<std::string>("Boundary condition types", 2);
	for(int i=0; i<2; i++) {
		if(types[i]=="Transparent")
			_boundary_types.push_back(BoundaryType::Transparent);
		else if (types[i]=="Periodic")
			_boundary_types.push_back(BoundaryType::Periodic);
		else
			throw std::string(
				"Parsing error: each item in \"Boundary condition types\" \n"
				"must be \"Transparent\" or \"Periodic\"");
	}
}


AlgorithmSettings::AlgorithmSettings(const nlohmann::json &j) :
	BaseSettings(j, "Algorithm settings"),
	general(json_node),
	bpm(json_node) {}
