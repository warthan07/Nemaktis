#ifndef MATERIALPROPERTIES_H
#define MATERIALPROPERTIES_H

#include "json.h"
using json = nlohmann::json;

class MaterialProperties {
public:
	MaterialProperties(json &material_settings);

	const double n_out, n_extra, n_ordi, n_lc, n_host;
	const std::vector<double> n_lower_layers, n_upper_layers;
	const std::vector<double> e_lower_layers, e_upper_layers;
};

#endif
