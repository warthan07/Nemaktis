#include <memory>
#include <boost/regex.hpp>

#include "BPMIteration.h"
#include "ParaxialPrimaryEvolutionOperator.h"
#include "PhaseEvolutionOperator.h"
#include "FresnelOperator.h"

BPMIteration::BPMIteration(
		const VectorField<double> &lc_sol,
		std::vector<VectorField<std::complex<double> > > (&bpm_sol)[2],
		const PhysicsCoefficients &coefs,
		const RootSettings &settings) :
	coefs(coefs),
	settings(settings),
	delta_x(lc_sol.mesh.delta_x),
	delta_y(lc_sol.mesh.delta_y),
	delta_z(lc_sol.mesh.delta_z),
	Nx(lc_sol.mesh.Nx),
	Ny(lc_sol.mesh.Ny),
	Nz(lc_sol.mesh.Nz),
	Nz_substeps(settings.algorithm.bpm.Nz_substeps),
	lc_sol(lc_sol),
	bpm_sol(bpm_sol),
	wavelengths(coefs.wavelengths()) {

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
}

void BPMIteration::update_optical_field() {

	TransverseOpticalField Eperp[2] = {{lc_sol.mesh}, {lc_sol.mesh}};

	std::string pol_strs[2] = {"X", "Y"};
	for(int wave_idx=0; wave_idx<wavelengths.size(); wave_idx++) {
		std::cout <<
			"[ wavelength: " << wavelengths[wave_idx] << "µm ]" << std::endl <<
			"\tInitializing evolution operators..." << std::endl;

		PermittivityTensorField eps(lc_sol, coefs, wavelengths[wave_idx]);
		ParaxialPrimaryEvolutionOperator primary_evolution_operator(
			eps, wavelengths[wave_idx], settings);
		PhaseEvolutionOperator secondary_evolution_operator(
			eps, wavelengths[wave_idx], settings);
		FresnelOperator fresnel_operator(eps, coefs.get_nin(wavelengths[wave_idx]));
		
		primary_evolution_operator.z_reinit();
		secondary_evolution_operator.z_reinit();

		// We initialize the values of the optical field on the entrance
		// plane
		std::cout <<
			"\tInitializing optical fields..." << std::endl;
		for(int pol_idx=0; pol_idx<2; pol_idx++) {
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
				Eperp[pol_idx](iperp,0) = beam_profile->get_Ex(x, y);
				Eperp[pol_idx](iperp,1) = beam_profile->get_Ey(x, y);
			}
			fresnel_operator.apply(Eperp[pol_idx]);

			#pragma omp parallel for
			for(int iperp=0; iperp<Nx*Ny; iperp++)
				for(int comp=0; comp<2; comp++)
					bpm_sol[pol_idx][wave_idx](iperp, comp) = Eperp[pol_idx](iperp,comp);
		}
		

		// We propagate the fields with the evolution operators.
		std::cout << "\tPropagating optical fields..." << std::endl;
		for(int iz=0; iz<Nz-1; iz++) {
			// We use the evolution operators to compute the optical field in
			// the next transverse plane.
			primary_evolution_operator.update();
			secondary_evolution_operator.update();

			for(int pol_idx=0; pol_idx<2; pol_idx++) {
				for(unsigned int k=0; k<Nz_substeps; k++) {
					secondary_evolution_operator.apply(Eperp[pol_idx]);
					primary_evolution_operator.apply(Eperp[pol_idx]);
					secondary_evolution_operator.apply(Eperp[pol_idx]);
				}
			}

			primary_evolution_operator.z_step_increment();
			secondary_evolution_operator.z_step_increment();
		
			// We save the calculated optical field values of the
			// new transverse plane in the solution vector
			int global_index;
			#pragma omp parallel for firstprivate(global_index)
			for(int iperp=0; iperp<Nx*Ny; iperp++) {
				global_index = (iperp+Nx*Ny*(iz+1));
				for(int pol_idx=0; pol_idx<2; pol_idx++) {
					for(int comp=0; comp<2; comp++)
						bpm_sol[pol_idx][wave_idx](global_index, comp) =
							Eperp[pol_idx](iperp,comp);
				}
			}
		}
	}
}
