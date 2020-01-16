#include "ScreenOutput.h"
#include "VolumeOutput.h"

#include "RootSettings.h"
#include "VTIReader.h"
#include "BPMIteration.h"

RootSettings::RootSettings(const nlohmann::json &j) :
	BaseSettings(j, ""),
	algorithm(j),
	physics(j),
	postprocessor(j) {}

void runFromSettings(RootSettings &settings) {

	/***********************************
	 * Solution vectors initialization *
	 ***********************************/

	// We initialize the director and bpm vectors
	std::cout <<
		std::endl << "Setting up the initial solution..." << std::endl;

	std::string initial_solution_file =
		settings.physics.initial_conditions.initial_solution_file;
	std::shared_ptr<VectorField<double> > lc_sol;

	if(initial_solution_file.substr(initial_solution_file.size()-4, 4) == ".vti") {
		VTIReader vti_reader(
			initial_solution_file,
			(settings.algorithm.general.lc_field_type=="Director") ? 3 : 6);
		vti_reader.fill_solution_vector(
			lc_sol, settings.physics.initial_conditions.basis_convention());
		lc_sol->set_mask_from_nonzeros();
	}
	else
		throw std::string(
			"Unrecognised initial solution filetype: must be a path to a vti file.");

	PhysicsCoefficients coefs(settings, lc_sol->mesh);
	std::vector<VectorField<std::complex<double> > >  bpm_sol[2] = {
		{coefs.wavelengths().size(), VectorField<std::complex<double> >(lc_sol->mesh, 2)},
		{coefs.wavelengths().size(), VectorField<std::complex<double> >(lc_sol->mesh, 2)}};


	/*********************************
	 * Postprocessors Initialization *
	 *********************************/

	typedef std::vector<std::shared_ptr<PostprocessorBase> > PostprocessorVec;
	
	PostprocessorVec postprocessors;
	PostprocessorSettings& pp_stg =
		settings.postprocessor;

	if(pp_stg.volume_output.activate)
		postprocessors.push_back(std::make_shared<VolumeOutput>(settings, coefs));
	if(pp_stg.micrograph_output.activate)
		postprocessors.push_back(std::make_shared<ScreenOutput>(settings, coefs));

	
	/*******************
	 * Problem solving * 
	 *******************/
	
	auto bpm_iteration = new BPMIteration(*lc_sol, bpm_sol, coefs, settings);

	bpm_iteration->update_optical_field();
	for(auto postprocessor : postprocessors) {
		postprocessor->apply(*lc_sol,bpm_sol);
	}

	postprocessors.clear();
}
