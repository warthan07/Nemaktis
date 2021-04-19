#include <string.h>
#include <stdlib.h>

#include <vtkXMLImageDataReader.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkStringArray.h>
#include <vtkImageData.h>
#include <vtkFieldData.h>

#include "PhysicsSettings.h"

#define PI 3.1415926535897932

FieldSetupSettings::FieldSetupSettings(
		const nlohmann::json &j) :
	BaseSettings(j, "Field setup"),
	beam_profile_type(parse<std::string>("Beam profile")),
	sample_file(parse<std::string>("Microscope sample file")) {

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

CoefficientsSettings::CoefficientsSettings(
		const nlohmann::json &j,
		const std::string &sample_file) :
	BaseSettings(j, "Coefficients"),
	nin_expression(
		(json_node.find("nin") != json_node.end()) ?
		parse<std::string>("nin") : "1"),
	nout_expression(
		(json_node.find("nout") != json_node.end()) ?
		parse<std::string>("nout") : "1") {
		
	if(sample_file!="") {
		auto reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
		reader->SetFileName(sample_file.c_str());
		reader->UpdateInformation();
		for(int i=0; i<reader->GetNumberOfPointArrays(); i++)
			reader->SetPointArrayStatus(reader->GetPointArrayName(i), 0);
		auto field_data = reader->GetOutput()->GetFieldData();

		if(field_data->GetAbstractArray("domain_refractive_indices")==nullptr)
			throw std::string(
				"Missing field array 'domain_refractive_indices' in vti sample file");
		auto nsample_data = vtkStringArray::SafeDownCast(
			field_data->GetAbstractArray("domain_refractive_indices"));
		for(int idx=0; idx<nsample_data->GetNumberOfValues(); idx++)
			_nsample_expressions.push_back(nsample_data->GetValue(idx));

		if(field_data->GetAbstractArray("up_iso_layer_refractive_indices")!=nullptr) {
			auto niso_up_data = vtkStringArray::SafeDownCast(
				field_data->GetAbstractArray("up_iso_layer_refractive_indices"));
			for(int idx=0; idx<niso_up_data->GetNumberOfValues(); idx++) 
				_niso_up_expressions.push_back(niso_up_data->GetValue(idx));
		}
		if(field_data->GetAbstractArray("up_iso_layer_thicknesses")!=nullptr) {
			auto hiso_up_data = vtkDoubleArray::FastDownCast(
				field_data->GetAbstractArray("up_iso_layer_thicknesses"));
			for(int idx=0; idx<hiso_up_data->GetNumberOfValues(); idx++) 
				_hiso_up.push_back(hiso_up_data->GetValue(idx));
		}
		if(field_data->GetAbstractArray("lo_iso_layer_refractive_indices")!=nullptr) {
			auto niso_lo_data = vtkStringArray::SafeDownCast(
				field_data->GetAbstractArray("lo_iso_layer_refractive_indices"));
			for(int idx=0; idx<niso_lo_data->GetNumberOfValues(); idx++) 
				_niso_lo_expressions.push_back(niso_lo_data->GetValue(idx));
		}
		if(field_data->GetAbstractArray("lo_iso_layer_thicknesses")!=nullptr) {
			auto hiso_lo_data = vtkDoubleArray::FastDownCast(
				field_data->GetAbstractArray("lo_iso_layer_thicknesses"));
			for(int idx=0; idx<hiso_lo_data->GetNumberOfValues(); idx++) 
				_hiso_lo.push_back(hiso_lo_data->GetValue(idx));
		}
	}
	else {
		_nsample_expressions = parse_vector<std::string>("nsample");
		if(json_node.find("niso_up") != json_node.end())
			_niso_up_expressions = parse_vector<std::string>("niso_up");
		if(json_node.find("niso_lo") != json_node.end())
			_niso_lo_expressions = parse_vector<std::string>("niso_lo");
		if(json_node.find("hiso_up") != json_node.end())
			_hiso_up = parse_vector<double>("hiso_up");
		if(json_node.find("hiso_lo") != json_node.end())
			_hiso_lo = parse_vector<double>("hiso_lo");
	}
	if(_nsample_expressions.size()%3 != 0)
		throw std::string(
			"Wrong number of values in array 'nsample'");
	if(_niso_up_expressions.size()!=_hiso_up.size())
		throw std::string(
			"Incorrect specification of indices/thicknesses of upper isotropic layers");
	if(_niso_lo_expressions.size()!=_hiso_lo.size())
		throw std::string(
			"Incorrect specification of indices/thicknesses of lower isotropic layers");

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
				double phi = iphi*PI/(3*ir);
				_q_vals.push_back(std::pair<double,double>(
					ir*NA/(Nc-1)*std::cos(phi),
					ir*NA/(Nc-1)*std::sin(phi)));
			}
		}
	}
}

PhysicsSettings::PhysicsSettings(const nlohmann::json &j) :
	BaseSettings(j, "Physics settings"),
	field_setup(json_node),
	coefs(json_node,field_setup.sample_file) {}
