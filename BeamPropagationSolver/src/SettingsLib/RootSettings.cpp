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

	std::cout <<
		std::endl << "Setting up the initial solution..." << std::endl;

	std::string initial_solution_file =
		settings.physics.initial_conditions.initial_solution_file;
	std::shared_ptr<VectorField<double> > lc_field;

	if(initial_solution_file.substr(initial_solution_file.size()-4, 4) == ".vti") {
		VTIReader vti_reader(
			initial_solution_file,
			(settings.algorithm.general.lc_field_type=="Director") ? 3 : 6);
		vti_reader.fill_solution_vector(
			lc_field, settings.physics.initial_conditions.basis_convention());
		lc_field->set_mask_from_nonzeros();
	}
	else
		throw std::string(
			"Unrecognised initial solution filetype: must be a path to a vti file.");

	PhysicsCoefficients coefs(settings, lc_field->mesh);
	ScreenOpticalFieldCollection screen_optical_fields(lc_field->mesh, coefs);


	/*******************
	 * Problem solving * 
	 *******************/
	
	BPMIteration bpm_iteration(*lc_field, screen_optical_fields, coefs, settings);
	bpm_iteration.propagate_fields();

	if(settings.postprocessor.volume_output.activate) {
		VolumeOutput volume_output(settings, coefs);
		volume_output.apply(*bpm_iteration.get_bulk_optical_fields());
	}
	if(settings.postprocessor.micrograph_output.activate) {
		ScreenOutput micrograph_output(settings, coefs);
		micrograph_output.apply(screen_optical_fields);
	}
}
