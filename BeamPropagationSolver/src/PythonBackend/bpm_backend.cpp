#include "BPMIteration.h"
#include "ScreenOutput.h"
#include "VolumeOutput.h"

void run_backend_without_mask(
		std::string json_str,
		double* director_vals, int n_director,
		std::complex<double>* fields_vals, int n_fields) {

	try {
		std::cout <<
			"[TESTING BRANCH]" << std::endl;
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
		auto nfield = std::make_shared<VectorField<double> >(
			mesh, 3, director_vals, n_director);

		PhysicsCoefficients coefs(settings, nfield->mesh);
		ScreenOpticalFieldCollection screen_optical_fields(nfield->mesh, coefs);

		BPMIteration bpm_iteration(
			*nfield, screen_optical_fields, coefs, settings);
		bpm_iteration.propagate_fields();

		ScreenOutput micrograph_output(settings, coefs);
		if(n_fields != coefs.wavelengths().size()*coefs.q_vals().size()*4*dims[0]*dims[1])
			throw std::string("Wrong dimension for the raw fields array");
		micrograph_output.apply_no_export(screen_optical_fields, fields_vals);

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
		std::string json_str,
		double* director_vals, int n_director,
		double* mask_vals, int n_mask,
		std::complex<double>* fields_vals, int n_fields) {

	try {
		std::cout <<
			"[TESTING BRANCH]" << std::endl;
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
		auto nfield = std::make_shared<VectorField<double> >(
			mesh, 3, director_vals, n_director, mask_vals, n_mask);

		PhysicsCoefficients coefs(settings, nfield->mesh);
		ScreenOpticalFieldCollection screen_optical_fields(nfield->mesh, coefs);

		BPMIteration bpm_iteration(
			*nfield, screen_optical_fields, coefs, settings);
		bpm_iteration.propagate_fields();

		ScreenOutput micrograph_output(settings, coefs);
		if(n_fields != coefs.wavelengths().size()*coefs.q_vals().size()*4*dims[0]*dims[1])
			throw std::string("Wrong dimension for the raw fields array");
		micrograph_output.apply_no_export(screen_optical_fields, fields_vals);

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
