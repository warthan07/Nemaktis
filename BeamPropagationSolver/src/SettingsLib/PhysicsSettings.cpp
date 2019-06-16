#include <string.h>
#include <stdlib.h>

#include "PhysicsSettings.h"

InitialConditionsSettings::InitialConditionsSettings(
		const nlohmann::json &j) :
	BaseSettings(j, "Initial conditions"),
	beam_profile_type(parse<std::string>("Beam profile")),
	initial_solution_file(parse<std::string>("LC field file")) {

	auto basis = parse<std::string>("Basis convention");
	if(basis=="XYZ")
		_basis_convention = BasisConvention::XYZ;
	else if(basis=="YZX")
		_basis_convention = BasisConvention::YZX;
	else if(basis=="ZXY")
		_basis_convention = BasisConvention::ZXY;
	else
		throw std::string(
			"Parsing error: \"Basis convention\" must be \"XYZ\", \"YZX\" or \"ZXY\"");
}	

CoefficientsSettings::CoefficientsSettings(const nlohmann::json &j) :
	BaseSettings(j, "Coefficients"),
	no_expression(parse<std::string>("no")),
	ne_expression(parse<std::string>("ne")),
	nhost_expression(
		(json_node.find("nhost") != json_node.end()) ?
		parse<std::string>("nhost") : "1") {

	if(json_node.find("Wavelengths") != json_node.end())
		_wavelengths = parse_vector<double>("Wavelengths");
	else {
		double mean_wavelength = parse<double>("Mean wavelength");
		double spectral_fwhm = parse<double>("Spectral FWHM");
		unsigned int N_wavelengths = parse<unsigned int>("N wavelengths");

		if(N_wavelengths==1)
			_wavelengths.push_back(mean_wavelength);
		else {
			_wavelengths.resize(N_wavelengths);
			for(unsigned int i=0; i<_wavelengths.size(); i++)
				_wavelengths[i] =
					mean_wavelength + spectral_fwhm*(double(i)/(_wavelengths.size()-1.) - 0.5);
		}
	}
}

PhysicsSettings::PhysicsSettings(const nlohmann::json &j) :
	BaseSettings(j, "Physics settings"),
	initial_conditions(json_node),
	coefs(json_node) {}
