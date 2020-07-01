#include <string.h>
#include <stdlib.h>

#include "PhysicsSettings.h"

#define PI 3.1415926535897932

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
		parse<std::string>("nhost") : "1"),
	nin_expression(
		(json_node.find("nin") != json_node.end()) ?
		parse<std::string>("nin") : "1"),
	nout_expression(
		(json_node.find("nout") != json_node.end()) ?
		parse<std::string>("nout") : "1") {

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

	if(json_node.find("Wavevectors") != json_node.end()) {
		auto raw_vals = parse_vector<double>("Wavevectors");
		if(raw_vals.size()%2!=0)
			throw std::string(
				"Wrong number of values for the \"Wavevectors\" entry,\n"
				"should be a multiple of 2");
		for(int q_idx=0; q_idx<(raw_vals.size()/2); q_idx++)
			_q_vals.push_back(std::pair<double,double>(
				raw_vals[2*q_idx], raw_vals[2*q_idx+1]));
	}
	else {
		double NA = parse<double>("Condenser numerical aperture");
		int Nc = parse<int>("N radial illumination directions");
		_q_vals.push_back(std::pair<double,double>(0,0));
		for(int ir=1; ir<Nc; ir++) {
			for(int iphi=0; iphi<6*ir; iphi++) {
				double phi = iphi*PI/3;
				_q_vals.push_back(std::pair<double,double>(
					ir*NA/(Nc-1)*std::cos(phi),
					ir*NA/(Nc-1)*std::sin(phi)));
			}
		}
	}
}

PhysicsSettings::PhysicsSettings(const nlohmann::json &j) :
	BaseSettings(j, "Physics settings"),
	initial_conditions(json_node),
	coefs(json_node) {}
