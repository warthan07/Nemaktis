#include <memory>
#include <boost/regex.hpp>

#include "BPMIteration.h"
#include "ParaxialPrimaryEvolutionOperator.h"
#include "PhaseEvolutionOperator.h"
#include "FresnelOperator.h"
#include "SimpleShiftOperator.h"

BPMIteration::BPMIteration(
		const VectorField<double> &lc_field,
		ScreenOpticalFieldCollection &screen_optical_fields,
		const PhysicsCoefficients &coefs,
		const RootSettings &settings) :
	coefs(coefs),
	settings(settings),
	delta_x(lc_field.mesh.delta_x),
	delta_y(lc_field.mesh.delta_y),
	delta_z(lc_field.mesh.delta_z),
	Nx(lc_field.mesh.Nx),
	Ny(lc_field.mesh.Ny),
	Nz(lc_field.mesh.Nz),
	Nz_substeps(settings.algorithm.bpm.Nz_substeps),
	lc_field(lc_field),
	screen_optical_fields(screen_optical_fields),
	bulk_output(settings.postprocessor.volume_output.activate),
	wavelengths(coefs.wavelengths()),
	q_vals(coefs.q_vals()) {

	auto& ic_settings = settings.physics.initial_conditions;
	if(ic_settings.beam_profile_type == "UniformBeam")
		beam_profile_type = BeamProfileType::UniformBeam;
	else {
		boost::smatch match_res;
		if(boost::regex_match(
				ic_settings.beam_profile_type, match_res,
				boost::regex("GaussianBeam\\(([0-9]+\\.?[0-9]*)\\)"))) {
			beam_profile_type = BeamProfileType::GaussianBeam;
			waist = std::stod(match_res[1]);
		}
		else {
			throw std::string("Unrecognised beam profile");
		}
	}

	if(bulk_output)
		bulk_optical_fields = std::make_shared<BulkOpticalFieldCollection>(
			lc_field.mesh, coefs);
}

void BPMIteration::propagate_fields() {

	for(int wave_idx=0; wave_idx<wavelengths.size(); wave_idx++) {
		std::cout <<
			"[ wavelength: " << wavelengths[wave_idx] << "µm ]" << std::endl <<
			"\tInitializing evolution operators..." << std::endl;

		PermittivityTensorField eps(lc_field, coefs, wavelengths[wave_idx]);
		ParaxialPrimaryEvolutionOperator primary_evolution_operator(
			eps, wavelengths[wave_idx], settings);
		PhaseEvolutionOperator secondary_evolution_operator(
			eps, wavelengths[wave_idx], settings);
		InputFresnelOperator input_fresnel_operator(eps, coefs.get_nin(wavelengths[wave_idx]));

		double nout_fresnel;
		auto& micrograph_output = settings.postprocessor.micrograph_output;
		if(micrograph_output.iso_layer_index.size()>0)
			nout_fresnel = micrograph_output.iso_layer_index[0];
		else
			nout_fresnel = coefs.get_nout(wavelengths[wave_idx]);
		OutputFresnelOperator output_fresnel_operator(eps, nout_fresnel);

		std::vector<SimpleShiftOperator> shift_operators;
		for(int q_idx=0; q_idx<q_vals.size(); q_idx++)
			shift_operators.emplace_back(eps, wavelengths[wave_idx], q_vals[q_idx], settings);
		
		primary_evolution_operator.z_reinit();
		secondary_evolution_operator.z_reinit();

		// We initialize the values of the optical field on the entrance
		// plane
		std::cout <<
			"\tInitializing optical fields..." << std::endl;
		for(int pol_idx=0; pol_idx<2; pol_idx++) {
			for(int q_idx=0; q_idx<q_vals.size(); q_idx++) {
				std::shared_ptr<BeamProfile> beam_profile;
				if(beam_profile_type == BeamProfileType::GaussianBeam)
					beam_profile = std::make_shared<GaussianBeam>(
						coefs, waist, 0.5*PI*pol_idx, wavelengths[wave_idx]);
				else if(beam_profile_type == BeamProfileType::UniformBeam)
					beam_profile = std::make_shared<UniformBeam>(
						coefs, 0.5*PI*pol_idx, wavelengths[wave_idx]);
	
				double x, y;
				for(int iperp=0; iperp<Nx*Ny; iperp++) {
					x = ((iperp%Nx) - (Nx-1)/2.) * delta_x;
					y = ((iperp/Nx) - (Ny-1)/2.) * delta_y;
					screen_optical_fields(wave_idx,q_idx,pol_idx)(iperp,0) =
						beam_profile->get_Ex(x, y);
					screen_optical_fields(wave_idx,q_idx,pol_idx)(iperp,1) =
						beam_profile->get_Ey(x, y);
				}
				input_fresnel_operator.apply(
					screen_optical_fields(wave_idx,q_idx,pol_idx));

				if(bulk_output) {
					#pragma omp parallel for
					for(int iperp=0; iperp<Nx*Ny; iperp++)
						for(int comp=0; comp<2; comp++)
							bulk_optical_fields->set_field_val(
								wave_idx, q_idx, pol_idx, comp, iperp,
								screen_optical_fields(wave_idx,q_idx,pol_idx)(iperp,comp));
				}
			}
		}
		

		// We propagate the fields with the evolution operators.
		std::cout << "\tPropagating optical fields..." << std::endl;
		for(int iz=0; iz<Nz-1; iz++) {
			// We use the evolution operators to compute the optical field in
			// the next transverse plane.
			primary_evolution_operator.update();
			secondary_evolution_operator.update();

			for(int pol_idx=0; pol_idx<2; pol_idx++) {
				for(int q_idx=0; q_idx<q_vals.size(); q_idx++) {
					shift_operators[q_idx].apply(
						screen_optical_fields(wave_idx,q_idx,pol_idx));
					secondary_evolution_operator.apply(
						screen_optical_fields(wave_idx,q_idx,pol_idx));
					primary_evolution_operator.apply(
						screen_optical_fields(wave_idx,q_idx,pol_idx));
					secondary_evolution_operator.apply(
						screen_optical_fields(wave_idx,q_idx,pol_idx));
					shift_operators[q_idx].apply(
						screen_optical_fields(wave_idx,q_idx,pol_idx));
				}
			}

			primary_evolution_operator.z_step_increment();
			secondary_evolution_operator.z_step_increment();
		
			// We save the calculated optical field values of the new transverse plane in
			// the bulk output object if the user asked for it.
			if(bulk_output) {
				int global_index;
				#pragma omp parallel for firstprivate(global_index)
				for(int iperp=0; iperp<Nx*Ny; iperp++) {
					global_index = (iperp+Nx*Ny*(iz+1));
					for(int pol_idx=0; pol_idx<2; pol_idx++) {
						for(int q_idx=0; q_idx<q_vals.size(); q_idx++) {
							for(int comp=0; comp<2; comp++)
								bulk_optical_fields->set_field_val(
									wave_idx, q_idx, pol_idx, comp, global_index,
									screen_optical_fields(wave_idx,q_idx,pol_idx)(iperp,comp));
						}
					}
				}
			}
		}

		// We apply the output Fresnel boundary conditions
		for(int pol_idx=0; pol_idx<2; pol_idx++)
			for(int q_idx=0; q_idx<q_vals.size(); q_idx++)
				output_fresnel_operator.apply(
					screen_optical_fields(wave_idx,q_idx,pol_idx));
	}
}
