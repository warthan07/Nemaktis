#include <memory>
#include <boost/regex.hpp>

#include "BPMIteration.h"
#include "ADIOperatorX.h"
#include "ADIOperatorY.h"
#include "PhaseOperator.h"
#include "FresnelOperator.h"

BPMIteration::BPMIteration(
		const VectorField<double> &lc_sol,
		std::vector<VectorField<std::complex<double> > > (&bpm_sol)[2],
		const PhysicsCoefficients &coefs,
		const RootSettings &settings) :
	coefs(coefs),
	bpm_settings(settings.algorithm.bpm),
	delta_x(lc_sol.mesh.delta_x),
	delta_y(lc_sol.mesh.delta_y),
	delta_z(lc_sol.mesh.delta_z),
	Nx(lc_sol.mesh.Nx),
	Ny(lc_sol.mesh.Ny),
	Nz(lc_sol.mesh.Nz),
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

	CartesianMesh transverse_mesh(
		{lc_sol.mesh.delta_x, lc_sol.mesh.delta_y, 0},
		{lc_sol.mesh.Nx, lc_sol.mesh.Ny, 1});
	VectorField<std::complex<double> > transverse_field_1(transverse_mesh, 2);
	VectorField<std::complex<double> > transverse_field_2(transverse_mesh, 2);

	std::string pol_strs[2] = {"X", "Y"};
	for(int wave_idx=0; wave_idx<wavelengths.size(); wave_idx++) {
		std::cout <<
			"[ wavelength: " << wavelengths[wave_idx] << "µm ]" << std::endl <<
			"\tInitializing permittivity and ADI/phase operators..." << std::endl;

		PermittivityTensorField eps(lc_sol, coefs, wavelengths[wave_idx]);
		ADIOperatorX adi_operator_x(eps, wavelengths[wave_idx], bpm_settings);
		ADIOperatorY adi_operator_y(eps, wavelengths[wave_idx], bpm_settings);
		PhaseOperator phase_operator(eps, wavelengths[wave_idx]);
		FresnelOperator fresnel_operator(eps, coefs.get_nin(wavelengths[wave_idx]));
		
		for(int pol_idx=0; pol_idx<2; pol_idx++) {
			adi_operator_x.z_step_reset();
			adi_operator_y.z_step_reset();
			phase_operator.z_step_reset();

			std::cout <<
				"\t[ input polarisation: " << pol_strs[pol_idx] << " ]" << std::endl;

			// We initialize the values of the optical field on the entrance
			// plane
			std::cout <<
				"\t\tInitializing optical fields..." << std::endl;
		
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
				transverse_field_2(iperp,0) = beam_profile->get_Ex(x, y);
				transverse_field_2(iperp,1) = beam_profile->get_Ey(x, y);
			}
			fresnel_operator.vmult(transverse_field_1, transverse_field_2);
		
			#pragma omp parallel for
			for(int iperp=0; iperp<Nx*Ny; iperp++)
				for(int comp=0; comp<2; comp++)
					bpm_sol[pol_idx][wave_idx](iperp, comp) = transverse_field_1(iperp,comp);

			// We propagate the fields with the ADI operators.
			std::cout << "\t\tPropagating optical fields..." << std::endl;
			for(int iz=0; iz<Nz-1; iz++) {
				// We use the phase and adi operators to compute the optical field in
				// the next transverse plane.
				adi_operator_x.update_tbc_wavectors(transverse_field_1);
				adi_operator_y.update_tbc_wavectors(transverse_field_1);

				phase_operator.vmult(transverse_field_2, transverse_field_1);

				adi_operator_y.switch_to_forward_operator();
				adi_operator_y.vmult(transverse_field_1, transverse_field_2);

				adi_operator_x.switch_to_backward_operator();
				adi_operator_x.inverse_vmult(transverse_field_2, transverse_field_1);

				adi_operator_x.switch_to_forward_operator();
				adi_operator_x.vmult(transverse_field_1, transverse_field_2);

				adi_operator_y.switch_to_backward_operator();
				adi_operator_y.inverse_vmult(transverse_field_2, transverse_field_1);

				phase_operator.vmult(transverse_field_1, transverse_field_2);

				adi_operator_x.z_step_increment();
				adi_operator_y.z_step_increment();
				phase_operator.z_step_increment();
		
				// We save the calculated optical field values of the
				// new transverse plane in the solution vector
				int global_index;
				#pragma omp parallel for firstprivate(global_index)
				for(int iperp=0; iperp<Nx*Ny; iperp++) {
					global_index = (iperp+Nx*Ny*(iz+1));
					for(int comp=0; comp<2; comp++)
						bpm_sol[pol_idx][wave_idx](global_index, comp) =
							transverse_field_1(iperp,comp);
				}
			}
		}
	}
}
