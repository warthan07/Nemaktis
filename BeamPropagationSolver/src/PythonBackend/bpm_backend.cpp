#include "BPMIteration.h"
#include "ScreenOutput.h"
#include "VolumeOutput.h"

void run_backend(
		std::string json_str,
		double* director_vals, int n_director,
		std::complex<double>* fields_vals, int n_fields) {

	try {
		auto j = nlohmann::json::parse(json_str);
		RootSettings settings(j);

		std::cout <<
			"Setting up the director field..." << std::endl;
	
		std::vector<double> spacings =
			j.at("Physics settings").at("Initial conditions").at("Mesh spacings");
		std::vector<int> dims =
			j.at("Physics settings").at("Initial conditions").at("Mesh dimensions");
		CartesianMesh mesh(
			{spacings[0], spacings[1], spacings[2]},
			{dims[0], dims[1], dims[2]});
		auto nfield = std::make_shared<VectorField<double> >(mesh, 3, director_vals, n_director);

		PhysicsCoefficients coefs(settings, nfield->mesh);
		std::vector<VectorField<std::complex<double> > >  bpm_sol[2] = {
			{coefs.wavelengths().size(), VectorField<std::complex<double> >(nfield->mesh, 2)},
			{coefs.wavelengths().size(), VectorField<std::complex<double> >(nfield->mesh, 2)}};

		auto bpm_iteration = new BPMIteration(*nfield, bpm_sol, coefs, settings);
		bpm_iteration->update_optical_field();

		ScreenOutput micrograph_output(settings, coefs);
		if(n_fields != coefs.wavelengths().size()*4*dims[0]*dims[1])
			throw std::string("Wrong dimension for the raw fields array");
		micrograph_output.apply_no_export(*nfield, bpm_sol, fields_vals);

		if(settings.postprocessor.volume_output.activate) {
			VolumeOutput volume_output(settings, coefs);
			volume_output.apply(*nfield, bpm_sol);
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
