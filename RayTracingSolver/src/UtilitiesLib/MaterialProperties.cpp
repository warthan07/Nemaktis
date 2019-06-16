#include "MaterialProperties.h"

MaterialProperties::MaterialProperties(json &j) :
	n_out(j.at("External medium refractive index")),
	n_extra(j.at("Extraordinary refractive index")),
	n_ordi(j.at("Ordinary refractive index")),
	n_lc((2*n_ordi+n_extra)/3.),
	n_host(j.at("Host fluid refractive index")),
	n_lower_layers(j.at("Refractive indices of the lower isotropic layers")
		.get<std::vector<double> >()),
	n_upper_layers(j.at("Refractive indices of the upper isotropic layers")
		.get<std::vector<double> >()),
	e_lower_layers(j.at("Thicknesses of the lower isotropic layers")
		.get<std::vector<double> >()),
	e_upper_layers(j.at("Thicknesses of the upper isotropic layers")
		.get<std::vector<double> >()) {

	if(n_lower_layers.size()!=e_lower_layers.size())
		throw(std::string(
			"Error: the number of lower isotropic layers must match when\n"
			"specifying the associated indices and thicknesses"));
	if(n_upper_layers.size()!=e_upper_layers.size())
		throw(std::string(
			"Error: the number of upper isotropic layers must match when\n"
			"specifying the associated indices and thicknesses"));
}
