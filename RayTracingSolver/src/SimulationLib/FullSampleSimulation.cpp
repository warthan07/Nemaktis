#include <fstream>
#include <sstream>

#include "FullSampleSimulation.h"
#include "LinearInterpolatedMapping.h"
#include "InverseScreenMap.h"

template <int dim>
FullSampleSimulation<dim>::FullSampleSimulation(
		json settings, double lc_thickness, unsigned int N_lc_steps,
		std::shared_ptr<CubicInterpolatedMapping<dim,3,double> > &n_field) :
	Simulation<dim>::Simulation(settings, lc_thickness, N_lc_steps, n_field) {
	
	auto N_rays_per_dim = parse_Vector<dim-1,int>(
		settings.at("Light source"), "Source N rays per dim");
	auto widths = parse_Vector<dim-1,double>(
		settings.at("Light source"), "Source widths");
	for(int i=0; i<dim-1; i++)
		widths(i) *= N_rays_per_dim(i)*1./(N_rays_per_dim(i)-4);
	auto origin = -widths/2;

	auto source_mesh =
		std::make_shared<CartesianMesh<dim-1> >(origin,widths,N_rays_per_dim);
	auto source_domain =
		std::make_shared<ParallelotopeDomain<dim-1> >(origin, widths);
	source = std::make_shared<LightSource<dim-1> >(
		settings, source_mesh, source_domain);

	extra_ray_family = std::make_shared<RayFamily<dim> >(
		*this->source, -this->lower_iso_thickness-this->lc_thickness, 0.5);
	ordi_ray_family = std::make_shared<RayFamily<dim> >(
		*this->source, -this->lower_iso_thickness-this->lc_thickness, 0.5);

	double eps_par = std::pow(this->mat_properties.n_extra, 2.);
	double eps_perp = std::pow(this->mat_properties.n_ordi, 2.);
	
	extra_ode = std::make_shared<ExtraordinaryODEFunction<dim> >(
		eps_par, eps_perp, n_field);
	ordi_ode = std::make_shared<OrdinaryODEFunction<dim> >(
		eps_perp, n_field);

	// We build the ODE sequences from bottom to top
	for(auto ode : this->lower_iso_odes) {
		extra_ray_ode_seq.odes.push_back(ode);
		ordi_ray_ode_seq.odes.push_back(ode);
	}

	extra_ray_ode_seq.odes.push_back(extra_ode);
	ordi_ray_ode_seq.odes.push_back(ordi_ode);
	
	for(auto ode : this->upper_iso_odes) {
		extra_ray_ode_seq.odes.push_back(ode);
		ordi_ray_ode_seq.odes.push_back(ode);
	}
	
	ray_splitting_ode_idx = int(this->lower_iso_odes.size()-1);

	N_rays = int(extra_ray_family->ray_bundles->size());
}

template <int dim>
void FullSampleSimulation<dim>::run() {

	std::vector<int> cur_ode_idx(N_rays, ray_splitting_ode_idx+1);

	auto e = *extra_ray_family;
	auto o = *ordi_ray_family;
	auto e_seq = extra_ray_ode_seq;
	auto o_seq = ordi_ray_ode_seq;

	RayStepper<dim> stepper(this->z_step);
	FresnelTransmission<dim> trans;
	std::shared_ptr<ODEFunction<dim> > trans_odes[4];
	std::shared_ptr<DefinitionDomain<dim-1> > screen_domain =
		std::make_shared<ParallelotopeDomain<dim-1> >(
			*(this->coarse_horiz_mesh));

	/****************************************************
	 * First, we propagate rays in the lower iso layers *
	 ****************************************************/
	std::cout <<
		"Initializing rays..." << std::endl << std::endl;
	ThreadException exc;

	std::vector<double> s_steps_extra(N_rays, this->z_step);
	std::vector<double> s_steps_ordi(N_rays, this->z_step);

	#pragma omp parallel for \
		firstprivate(e_seq, o_seq, stepper, trans) private(trans_odes)
	for(int i=0; i<N_rays; i++) {
		try {
			auto& extra_ray_bundle = e.ray_bundles->at(i);
			auto& ordi_ray_bundle = o.ray_bundles->at(i);

			for(int ode_idx=0; ode_idx<=ray_splitting_ode_idx; ode_idx++) {

				stepper.euler_step_to_boundary(
					e_seq.odes[ode_idx], extra_ray_bundle);
				stepper.euler_step_to_boundary(
					o_seq.odes[ode_idx], ordi_ray_bundle);

				trans_odes[0] = e_seq.odes[ode_idx];
				trans_odes[1] = o_seq.odes[ode_idx];
				trans_odes[2] = e_seq.odes[ode_idx+1];
				trans_odes[3] = o_seq.odes[ode_idx+1];
				
				// Special method for the ray splitting at the Iso/Aniso
				// interface
				if(ode_idx==ray_splitting_ode_idx)
					trans.update(trans_odes, extra_ray_bundle, ordi_ray_bundle);
				else {
					trans.update(trans_odes, extra_ray_bundle);
					trans.update(trans_odes, ordi_ray_bundle);
				}
			}
		}
		catch(...) {
			exc.capture_exception();
		}
	}
	exc.rethrow();

	/*************************************************
	 * Second, we propagate rays inside the LC layer *
	 *************************************************/
	double z;
	Vector<3,std::complex<double> > E_tot;
	std::complex<double> I(0., 1.);

	unsigned int node_idx = 0;
	for(int k=0; k<this->N_lc_steps; k++) {
		std::cout <<
			"Step " << k+1 << "/" << this->N_lc_steps << ":" << std::endl;

		if(this->bulk_fields_output || this->bulk_ray_output) {
			std::cout <<
				"	Updating map data..." << std::endl;
			e.update_map_data();
			o.update_map_data();
		}

		if(this->bulk_fields_output) {
			// We inverse the screen mapping
			std::cout <<
				"	Inversing ray mapping..." << std::endl;

			std::cout <<
				"		Coarse search of inverses..." << std::endl;

			InverseScreenMap<dim-1> inv_map(
				e.target_map, *this->coarse_horiz_mesh, screen_domain,
				this->typical_length, this->hc_parameters);
			for(int r=0; r<this->N_refinement_cycles; r++) {
				std::cout <<
					"		Refinement of inverses..." << std::endl;
				inv_map.refine_data();
			}
			auto inv_map_data = inv_map.get_data();


			// We reconstruct the field values and save them
			std::cout <<
				"	Reconstructing optical fields..." << std::endl;

			#pragma omp parallel for firstprivate(E_tot, e, o)
			for(int hpx_idx=0; hpx_idx<inv_map_data->size(); hpx_idx++) {
				auto it = inv_map_data->at(hpx_idx);
				int px_idx = hpx_idx + k*this->N_horiz_pixels;

				std::vector<Vector<3,std::complex<double> > > E_vals[2];
				std::vector<double> optical_length_vals;

				try {
					e.get_interpolated_maxwell_fields(
						it.second, E_vals, optical_length_vals);
					o.get_interpolated_maxwell_fields(
						it.first, E_vals, optical_length_vals);

					for(int wave_idx=0; wave_idx<this->N_wavelengths; wave_idx++) {
						for(int pol_idx=0; pol_idx<2; pol_idx++) {
							E_tot = std::complex<double>(0,0);
							for(int ray=0; ray<E_vals[0].size(); ray++)
								E_tot += E_vals[pol_idx][ray]*std::exp(
									I*2.*PI*(optical_length_vals[ray]
										-this->mat_properties.n_ordi*this->z_step*k) /
									this->wavelengths[wave_idx]);

							for(int comp=0; comp<3; comp++) {
								this->bulk_fields_data->E_real[pol_idx][wave_idx]
									->SetComponent(px_idx, comp, std::real(E_tot(comp)));
								this->bulk_fields_data->E_imag[pol_idx][wave_idx]
									->SetComponent(px_idx, comp, std::imag(E_tot(comp)));
							}
						}
					}
					this->bulk_fields_data->ray_multiplicity->SetComponent(
						px_idx, 0, E_vals[0].size());
				}
				catch(...) {
					exc.capture_exception();
				}
			}
			exc.rethrow();
		}

		if(this->bulk_ray_output) {
			std::cout <<
				"	Interpolating ray data..." << std::endl;

			#pragma omp parallel for firstprivate(E_tot, e, o)
			for(int hpx_idx=0; hpx_idx<this->N_horiz_pixels; hpx_idx++) {
				MultiDimIndex<dim-1> mesh_idx(
					this->full_horiz_mesh->n_points_per_dim, 0, false);
				mesh_idx = hpx_idx;
				Vector<dim-1,double> screen_pos = 
					this->full_horiz_mesh->origin +
					mesh_idx.get() * this->full_horiz_mesh->cell_lengths;

				int px_idx = hpx_idx + k*this->N_horiz_pixels;
				e.save_interpolated_ray_values(px_idx, screen_pos, this->bulk_extra_data);
				o.save_interpolated_ray_values(px_idx, screen_pos, this->bulk_ordi_data);
			}
		}


		// We evolve the rays to the next z-slice
		std::cout <<
			"	Propagating rays to the next target plane..." << std::endl << std::endl;

		if(k<this->N_lc_steps-1) {
			#pragma omp parallel for \
				firstprivate(e_seq,o_seq,stepper,trans,z) private(trans_odes)
			for(int i=0; i<N_rays; i++) {
				auto& extra_ray_bundle = e.ray_bundles->at(i);
				auto& ordi_ray_bundle = o.ray_bundles->at(i);
				int& ode_idx = cur_ode_idx[i];

				try {
					z = - this->lc_thickness/2 + (k+1)*this->z_step;

					stepper.euler_z_step(
						e_seq.odes[ode_idx], extra_ray_bundle, s_steps_extra[i], z);
					stepper.euler_z_step(
						o_seq.odes[ode_idx], ordi_ray_bundle, s_steps_ordi[i], z);
				}
				catch(...) {
					exc.capture_exception();
				}
			}
			exc.rethrow();
		}
		else {
			int ode_idx = ray_splitting_ode_idx+1;
			#pragma omp parallel for \
				firstprivate(e_seq, o_seq, stepper, trans) private(trans_odes)
			for(int i=0; i<N_rays; i++) {
				try {
					auto& extra_ray_bundle = e.ray_bundles->at(i);
					auto& ordi_ray_bundle = o.ray_bundles->at(i);

					stepper.euler_step_to_boundary(
						e_seq.odes[ode_idx], extra_ray_bundle);
					stepper.euler_step_to_boundary(
						o_seq.odes[ode_idx], ordi_ray_bundle);

					trans_odes[0] = e_seq.odes[ode_idx];
					trans_odes[1] = o_seq.odes[ode_idx];
					trans_odes[2] = e_seq.odes[ode_idx+1];
					trans_odes[3] = o_seq.odes[ode_idx+1];
					
					trans.update(trans_odes, extra_ray_bundle);
					trans.update(trans_odes, ordi_ray_bundle);
				}
				catch(...) {
					exc.capture_exception();
				}
			}
			exc.rethrow();
		}
	}
	if(this->bulk_ray_output) {
		std::cout <<
			"Exporting bulk ray data..." << std::endl;
		this->bulk_extra_data->write(
			this->basedir, this->bulk_basename+"_extra_rays.vti");
		this->bulk_ordi_data->write(
			this->basedir, this->bulk_basename+"_ordi_rays.vti");
	}
	if(this->bulk_fields_output) {
		std::cout <<
			"Exporting bulk optical fields..." << std::endl;
		this->bulk_fields_data->write(
			this->basedir, this->bulk_basename+"_fields.vti");
	}
	std::cout << std::endl;


	/***********************************************************
	 * Finally, we propagate rays through the upper iso layers *
	 ***********************************************************/

	// If we need reconstructed fields on the screen, we reconstruct
	// the fields locally on a centered focal plane before propagating
	// fields in Fourier space throught the iso layers
	if(this->screen_fields_output) {
		std::cout <<
			"Reconstructing optical fields on the focalisation plane..." <<
			std::endl;

		// We go to a centered virtual focal plane
		#pragma omp parallel for firstprivate(stepper)
		for(int i=0; i<N_rays; i++) {
			stepper.virtual_z_step(e.ray_bundles->at(i), this->z_foc_1);
			stepper.virtual_z_step(o.ray_bundles->at(i), this->z_foc_1);
		}

		std::cout <<
			"	Updating map data..." << std::endl;
		e.update_map_data();
		o.update_map_data();

		std::cout <<
			"	Inversing ray mapping..." << std::endl;

		std::cout <<
			"		Coarse search of inverses..." << std::endl;

		InverseScreenMap<dim-1> inv_map(
			e.target_map, *this->coarse_horiz_mesh, screen_domain,
			this->typical_length, this->hc_parameters);
		for(int r=0; r<this->N_refinement_cycles; r++) {
			std::cout <<
				"		Refinement of inverses..." << std::endl;
			inv_map.refine_data();
		}
		auto inv_map_data = inv_map.get_data();


		// We reconstruct the field values and save them
		std::cout <<
			"	Reconstructing optical fields locally..." << std::endl;

		#pragma omp parallel for firstprivate(E_tot, e, o)
		for(int hpx_idx=0; hpx_idx<inv_map_data->size(); hpx_idx++) {
			auto it = inv_map_data->at(hpx_idx);

			std::vector<Vector<3,std::complex<double> > > E_vals[2];
			std::vector<double> optical_length_vals;

			try {
				e.get_interpolated_maxwell_fields(
					it.second, E_vals, optical_length_vals);
				o.get_interpolated_maxwell_fields(
					it.first, E_vals, optical_length_vals);

				for(int wave_idx=0; wave_idx<this->N_wavelengths; wave_idx++) {
					for(int pol_idx=0; pol_idx<2; pol_idx++) {
						E_tot = std::complex<double>(0,0);
						for(int ray=0; ray<E_vals[0].size(); ray++)
							E_tot += E_vals[pol_idx][ray]*std::exp(
								I*2.*PI*optical_length_vals[ray] /
								this->wavelengths[wave_idx]);

						for(int comp=0; comp<2; comp++) {
							this->screen_fields_data->E_real[pol_idx][wave_idx]->SetComponent(
								hpx_idx, comp, std::real(E_tot(comp)));
							this->screen_fields_data->E_imag[pol_idx][wave_idx]->SetComponent(
								hpx_idx, comp, std::imag(E_tot(comp)));
						}
						this->screen_fields_data->E_real[pol_idx][wave_idx]->SetComponent(
							hpx_idx, 2, 0);
						this->screen_fields_data->E_imag[pol_idx][wave_idx]->SetComponent(
							hpx_idx, 2, 0);
					}
				}
				this->screen_fields_data->ray_multiplicity->SetComponent(
					hpx_idx, 0, E_vals[0].size());
			}
			catch(...) {
				exc.capture_exception();
			}
		}
		exc.rethrow();

		std::cout <<
			"	Switching to Fourier space to propagate through the iso layers..." <<
			std::endl;
		this->apply_fourier_iso_layer_filter();

		std::cout <<
			"	Exporting screen optical fields..." << std::endl << std::endl;
		this->screen_fields_data->write(
			this->basedir, this->screen_basename+"_fields.vti");

		// We come back inside the first iso layer
		#pragma omp parallel for firstprivate(stepper)
		for(int i=0; i<N_rays; i++) {
			stepper.virtual_z_step(e.ray_bundles->at(i), this->lc_thickness/2+1e-8);
			stepper.virtual_z_step(o.ray_bundles->at(i), this->lc_thickness/2+1e-8);
		}
	}
	

	// If we need ray data on the screen, we propagate rays in the usual
	// way throught the iso layers.
	if(this->screen_ray_output) {
		std::cout <<
			"Interpolating ray data on the focalisation plane..." << std::endl;
		#pragma omp parallel for \
		firstprivate(e_seq, o_seq, stepper, trans) private(trans_odes)
		for(int i=0; i<N_rays; i++) {
			try {
				auto& extra_ray_bundle = e.ray_bundles->at(i);
				auto& ordi_ray_bundle = o.ray_bundles->at(i);
	
				for(int ode_idx=ray_splitting_ode_idx+2; ode_idx<e_seq.odes.size()-1;
						ode_idx++) {
	
					stepper.euler_step_to_boundary(
						e_seq.odes[ode_idx], extra_ray_bundle);
					stepper.euler_step_to_boundary(
						o_seq.odes[ode_idx], ordi_ray_bundle);
	
					trans_odes[0] = e_seq.odes[ode_idx];
					trans_odes[1] = o_seq.odes[ode_idx];
					trans_odes[2] = e_seq.odes[ode_idx+1];
					trans_odes[3] = o_seq.odes[ode_idx+1];

					trans.update(trans_odes, extra_ray_bundle);
					trans.update(trans_odes, ordi_ray_bundle);
				}
			}
			catch(...) {
				exc.capture_exception();
			}
		}
		exc.rethrow();

		#pragma omp parallel for firstprivate(stepper)
		for(int i=0; i<N_rays; i++) {
			stepper.virtual_z_step(e.ray_bundles->at(i), this->z_foc_2);
			stepper.virtual_z_step(o.ray_bundles->at(i), this->z_foc_2);
		}

		std::cout <<
			"	Updating map data..." << std::endl;
		e.update_map_data();
		o.update_map_data();

		#pragma omp parallel for firstprivate(E_tot, e, o)
		for(int hpx_idx=0; hpx_idx<this->N_horiz_pixels; hpx_idx++) {
			MultiDimIndex<dim-1> mesh_idx(
				this->full_horiz_mesh->n_points_per_dim, 0, false);
			mesh_idx = hpx_idx;
			Vector<dim-1,double> screen_pos = 
				this->full_horiz_mesh->origin +
				mesh_idx.get() * this->full_horiz_mesh->cell_lengths;

			e.save_interpolated_ray_values(hpx_idx, screen_pos, this->screen_extra_data);
			o.save_interpolated_ray_values(hpx_idx, screen_pos, this->screen_ordi_data);
		}

		std::cout <<
			"	Exporting screen ray data..." << std::endl << std::endl;
		this->screen_extra_data->write(
			this->basedir, this->screen_basename+"_extra_rays.vti");
		this->screen_ordi_data->write(
			this->basedir, this->screen_basename+"_ordi_rays.vti");
	}
}

template class FullSampleSimulation<2>;
template class FullSampleSimulation<3>;
