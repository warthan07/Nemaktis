#include "BPMIteration.h"
#include "ScreenOutput.h"
#include "VolumeOutput.h"

void run_backend_without_mask(
		std::string json_str,
		double* lc_field_vals, int n_lc_vals,
		std::complex<double>* E_field_vals, int n_E_vals) {

	try {
		std::cout <<
			"Setting up the director field..." << std::endl;
	
		auto j = nlohmann::json::parse(json_str);
		RootSettings settings(j);

		std::vector<double> spacings =
			j.at("Physics settings").at("Initial conditions").at("Mesh spacings");
		std::vector<int> dims =
			j.at("Physics settings").at("Initial conditions").at("Mesh dimensions");
		CartesianMesh mesh(
			{spacings[0], spacings[1], spacings[2]},
			{dims[0], dims[1], dims[2]});
		auto lc_field = std::make_shared<VectorField<double> >(
			mesh, lc_field_vals, n_lc_vals);
		if(lc_field->field_dim!=3 && lc_field->field_dim!=6)
			throw std::string(
				"Unexpected dimension for the LC orientational field, should be 3 or 6.");

		PhysicsCoefficients coefs(settings, lc_field->mesh);
		ScreenOpticalFieldCollection screen_optical_fields(lc_field->mesh, coefs, E_field_vals, n_E_vals);

		BPMIteration bpm_iteration(
			*lc_field, screen_optical_fields, coefs, settings);
		bpm_iteration.propagate_fields();

		ScreenOutput micrograph_output(settings, coefs);
		micrograph_output.apply_no_export(screen_optical_fields);

		if(settings.postprocessor.volume_output.activate) {
			VolumeOutput volume_output(settings, coefs);
			volume_output.apply(*bpm_iteration.get_bulk_optical_fields());
		}
	}
	catch (std::exception &exc) {
		std::cerr <<
			std::endl << std::endl <<
			"-------------------------------------------------" <<
			std::endl <<
			"Exception on processing: " << std::endl <<
			exc.what() << std::endl <<
			"Aborting!" << std::endl <<
			"-------------------------------------------------" <<
			std::endl;
		return;
	}
	catch(std::string &str) {
		std::cerr <<
			std::endl << std::endl <<
			"-------------------------------------------------" <<
			std::endl <<
			"Exception on processing: " << std::endl <<
			str << std::endl <<
			"Aborting!" << std::endl <<
			"-------------------------------------------------" <<
			std::endl;
		return;
	}
	catch (...) {
		std::cerr <<
			std::endl << std::endl <<
			"-------------------------------------------------" <<
			std::endl <<
			"Unknown exception!" << std::endl <<
			"Aborting!" << std::endl <<
			"-------------------------------------------------" <<
			std::endl;
		return;
	}
}

void run_backend_with_mask(
		std::string json_str, std::string mask_formula,
		double* lc_field_vals, int n_lc_vals,
		double* mask_vals, int n_mask,
		std::complex<double>* E_field_vals, int n_E_vals) {

	try {
		std::cout <<
			"Setting up the director field..." << std::endl;
	
		auto j = nlohmann::json::parse(json_str);
		RootSettings settings(j);

		std::vector<double> spacings =
			j.at("Physics settings").at("Initial conditions").at("Mesh spacings");
		std::vector<int> dims =
			j.at("Physics settings").at("Initial conditions").at("Mesh dimensions");
		CartesianMesh mesh(
			{spacings[0], spacings[1], spacings[2]},
			{dims[0], dims[1], dims[2]});
		auto lc_field = std::make_shared<VectorField<double> >(
			mesh, lc_field_vals, n_lc_vals, mask_vals, n_mask, mask_formula);
		if(lc_field->field_dim!=3 && lc_field->field_dim!=6)
			throw std::string(
				"Unexpected dimension for the LC orientational field, should be 3 or 6.");

		PhysicsCoefficients coefs(settings, lc_field->mesh);
		ScreenOpticalFieldCollection screen_optical_fields(lc_field->mesh, coefs, E_field_vals, n_E_vals);

		BPMIteration bpm_iteration(
			*lc_field, screen_optical_fields, coefs, settings);
		bpm_iteration.propagate_fields();

		ScreenOutput micrograph_output(settings, coefs);
		if(n_E_vals != coefs.wavelengths().size()*coefs.q_vals().size()*4*dims[0]*dims[1])
			throw std::string("Wrong dimension for the raw fields array");
		micrograph_output.apply_no_export(screen_optical_fields);

		if(settings.postprocessor.volume_output.activate) {
			VolumeOutput volume_output(settings, coefs);
			volume_output.apply(*bpm_iteration.get_bulk_optical_fields());
		}
	}
	catch (std::exception &exc) {
		std::cerr <<
			std::endl << std::endl <<
			"-------------------------------------------------" <<
			std::endl <<
			"Exception on processing: " << std::endl <<
			exc.what() << std::endl <<
			"Aborting!" << std::endl <<
			"-------------------------------------------------" <<
			std::endl;
		return;
	}
	catch(std::string &str) {
		std::cerr <<
			std::endl << std::endl <<
			"-------------------------------------------------" <<
			std::endl <<
			"Exception on processing: " << std::endl <<
			str << std::endl <<
			"Aborting!" << std::endl <<
			"-------------------------------------------------" <<
			std::endl;
		return;
	}
	catch (...) {
		std::cerr <<
			std::endl << std::endl <<
			"-------------------------------------------------" <<
			std::endl <<
			"Unknown exception!" << std::endl <<
			"Aborting!" << std::endl <<
			"-------------------------------------------------" <<
			std::endl;
		return;
	}
}
