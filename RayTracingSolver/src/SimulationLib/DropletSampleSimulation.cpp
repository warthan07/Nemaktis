#include <fstream>
#include <sstream>

#include "DropletSampleSimulation.h"
#include "InverseScreenMap.h"

DropletSampleSimulation::DropletSampleSimulation(
		json settings, double lc_thickness, unsigned int N_lc_steps,
		std::shared_ptr<CubicInterpolatedMapping<3,3,double> > &n_field) :
	Simulation<3>::Simulation(settings, lc_thickness, N_lc_steps, n_field) {

	auto j = settings.at("Geometry").at("Droplet sample parameters");
	double dz = j.at("Distance from upper sample plate");
	if(dz<=0)
		throw(std::string(
			"The droplet is overlapping on the sample plates!"
			"Please specify a distance from the upper plate > 0"));
	droplet_radius = lc_thickness/2. - dz;

	auto N_rays_per_dim = parse_Vector<2,int>(
		settings.at("Light source"), "Source N rays per dim");
	auto widths = parse_Vector<2,double>(
		settings.at("Light source"), "Source widths");
	for(int i=0; i<2; i++)
		widths(i) *= N_rays_per_dim(i)*1./(N_rays_per_dim(i)-4);
	auto origin = -widths/2;

	auto source_mesh =
		std::make_shared<CartesianMesh<2> >(origin,widths,N_rays_per_dim);
	auto source_domain = std::make_shared<ParallelotopeDomain<2> >(origin, widths);

	// We restrict the radius for extra and ordi rays to avoid
	// numerical instability at the droplet boundary.
	double shrink_factor = j.at("Source shrink factor");
	if(shrink_factor<=0 || shrink_factor>=1)
		throw(std::string(
			"The shrink factor for the droplet source radius "
			"must be between 0 and 1"));
	double restricted_radius = shrink_factor * droplet_radius;

	auto restricted_droplet_source_domain = std::make_shared<SphericalDomain<2> >(
		Vector<2,double>(0), restricted_radius);

	auto droplet_source_domain = std::make_shared<SphericalDomain<2> >(
		Vector<2,double>(0), droplet_radius);
	auto iso_source_domain = std::make_shared<SubstractedDomain<2> >(
		source_domain, droplet_source_domain);

	iso_source = std::make_shared<LightSource<2> >(
		settings, source_mesh, iso_source_domain);
	droplet_source = std::make_shared<LightSource<2> >(
		settings, source_mesh, restricted_droplet_source_domain);

	extra_ray_family = std::make_shared<RayFamily<3> >(
		*droplet_source, -this->lower_iso_thickness-this->lc_thickness, 0.5);
	ordi_ray_family = std::make_shared<RayFamily<3> >(
		*droplet_source, -this->lower_iso_thickness-this->lc_thickness, 0.5);
	iso_ray_family = std::make_shared<RayFamily<3> >(
		*iso_source, -this->lower_iso_thickness-this->lc_thickness, 1);

	double eps_par = std::pow(this->mat_properties.n_extra, 2.);
	double eps_perp = std::pow(this->mat_properties.n_ordi, 2.);
	double eps_host = std::pow(this->mat_properties.n_host, 2.);

	auto sample_domain = std::make_shared<SlabDomain<3> >(
		BasisVector<3>(2), -0.5*lc_thickness, 0.5*lc_thickness);
	iso_ode = std::make_shared<IsotropicODEFunction<3> >(
		eps_host, std::make_shared<SubstractedDomain<3> >(
			sample_domain, n_field->get_def_domain()));
	extra_ode = std::make_shared<ExtraordinaryODEFunction<3> >(
		eps_par, eps_perp, n_field);
	ordi_ode = std::make_shared<OrdinaryODEFunction<3> >(
		eps_perp, n_field);

	// We build the ODE sequences from bottom to top
	for(auto ode : this->lower_iso_odes) {
		extra_ray_ode_seq.odes.push_back(ode);
		ordi_ray_ode_seq.odes.push_back(ode);
		iso_ray_ode_seq.odes.push_back(ode);
	}

	extra_ray_ode_seq.odes.push_back(iso_ode);
	extra_ray_ode_seq.odes.push_back(extra_ode);
	extra_ray_ode_seq.odes.push_back(iso_ode);

	ordi_ray_ode_seq.odes.push_back(iso_ode);
	ordi_ray_ode_seq.odes.push_back(ordi_ode);
	ordi_ray_ode_seq.odes.push_back(iso_ode);

	iso_ray_ode_seq.odes.push_back(iso_ode);

	for(auto ode : this->upper_iso_odes) {
		extra_ray_ode_seq.odes.push_back(ode);
		ordi_ray_ode_seq.odes.push_back(ode);
		iso_ray_ode_seq.odes.push_back(ode);
	}

	ray_splitting_ode_idx = int(this->lower_iso_odes.size());

	N_droplet_rays = int(extra_ray_family->ray_bundles->size());
	N_iso_rays = int(iso_ray_family->ray_bundles->size());
}

void DropletSampleSimulation::run() {

	std::vector<int> cur_e_ode_idx(N_droplet_rays, ray_splitting_ode_idx);
	std::vector<int> cur_o_ode_idx(N_droplet_rays, ray_splitting_ode_idx);

	auto e = *extra_ray_family;
	auto o = *ordi_ray_family;
	auto s = *iso_ray_family;

	auto e_seq = extra_ray_ode_seq;
	auto o_seq = ordi_ray_ode_seq;
	auto s_seq = iso_ray_ode_seq;

	RayStepper<3> stepper(this->z_step);
	FresnelTransmission<3> trans;
	std::shared_ptr<ODEFunction<3> > trans_odes[4];
	std::shared_ptr<DefinitionDomain<2> > screen_domain =
		std::make_shared<ParallelotopeDomain<2> >(
			*(this->coarse_horiz_mesh));


	/****************************************************
	 * First, we propagate rays in the lower iso layers *
	 ****************************************************/

	std::cout <<
		"Initializing rays..." << std::endl << std::endl;
	ThreadException exc;

	std::vector<double> s_steps_extra(N_droplet_rays, this->z_step);
	std::vector<double> s_steps_ordi(N_droplet_rays, this->z_step);
	std::vector<double> s_steps_iso(N_iso_rays, this->z_step);

	#pragma omp parallel for \
		firstprivate(e_seq, o_seq, stepper, trans) private(trans_odes)
	for(int ray_idx=0; ray_idx<N_droplet_rays; ray_idx++) {
		/*try {*/
			/*auto& extra_ray_bundle = e.ray_bundles->at(ray_idx);
			auto& ordi_ray_bundle = o.ray_bundles->at(ray_idx);

			for(int ode_idx=0; ode_idx<=ray_splitting_ode_idx-1; ode_idx++) {

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

			stepper.euler_z_step(
				e_seq.odes[ray_splitting_ode_idx], extra_ray_bundle,
				s_steps_extra[ray_idx], -droplet_radius);
			stepper.euler_z_step(
				o_seq.odes[ray_splitting_ode_idx], ordi_ray_bundle,
				s_steps_ordi[ray_idx], -droplet_radius);*/
		/*}
		catch(...) {
			exc.capture_exception();
		}*/
	}
	exc.rethrow();

	#pragma omp parallel for \
		firstprivate(s_seq, stepper, trans) private(trans_odes)
	for(int ray_idx=0; ray_idx<N_iso_rays; ray_idx++) {
		try {
			auto& iso_ray_bundle = s.ray_bundles->at(ray_idx);

			for(int ode_idx=0; ode_idx<=ray_splitting_ode_idx-1; ode_idx++) {

				stepper.euler_step_to_boundary(
					s_seq.odes[ode_idx], iso_ray_bundle);

				trans_odes[0] = s_seq.odes[ode_idx];
				trans_odes[1] = s_seq.odes[ode_idx];
				trans_odes[2] = s_seq.odes[ode_idx+1];
				trans_odes[3] = s_seq.odes[ode_idx+1];

				trans.update(trans_odes, iso_ray_bundle);
			}

			stepper.euler_z_step(
				s_seq.odes[ray_splitting_ode_idx], iso_ray_bundle,
				s_steps_iso[ray_idx], -droplet_radius);
		}
		catch(...) {
			exc.capture_exception();
		}
	}
	exc.rethrow();


	/*************************************************
	 * Second, we propagate rays inside the LC layer *
	 *************************************************/

	Vector<3,std::complex<double> > E_tot;
	std::complex<double> I(0., 1.);

	unsigned int node_idx = 0;
	for(int k=0; k<this->N_lc_steps; k++) {
		std::cout <<
			"Step " << k+1 << "/" << this->N_lc_steps << ":" << std::endl;

		if(this->bulk_fields_output || this->bulk_ray_output) {
			// We update the rays data
			std::cout <<
				"	Updating map data..." << std::endl;

			e.update_map_data();
			o.update_map_data();
			s.update_map_data();
		}

		if(this->bulk_fields_output) {
			// We inverse the screen mapping
			std::cout <<
				"	Inversing screen mapping..." << std::endl;

			std::cout <<
				"		Coarse search of inverses..." << std::endl;

			InverseScreenMap<2> e_inv_map(
				e.target_map, *this->coarse_horiz_mesh, screen_domain,
				this->typical_length, this->hc_parameters);
			InverseScreenMap<2> o_inv_map(
				o.target_map, *this->coarse_horiz_mesh, screen_domain,
				this->typical_length, this->hc_parameters);
			for(int r=0; r<this->N_refinement_cycles; r++) {
				std::cout <<
					"		Refinement of inverses..." << std::endl;
				e_inv_map.refine_data();
				o_inv_map.refine_data();
			}
			auto e_inv_map_data = e_inv_map.get_data();
			auto o_inv_map_data = o_inv_map.get_data();

			// We reconstruct the field values and save them
			std::cout <<
				"	Reconstructing optical fields..." << std::endl;

			#pragma omp parallel for \
				firstprivate(E_tot, e, o, s)
			for(int hpx_idx=0; hpx_idx<this->N_horiz_pixels; hpx_idx++) {
				auto e_it = e_inv_map_data->at(hpx_idx);
				auto o_it = o_inv_map_data->at(hpx_idx);
				int px_idx = hpx_idx + k*this->N_horiz_pixels;

				std::vector<Vector<3,std::complex<double> > > E_vals[2];
				std::vector<double> optical_length_vals;

				try {
					e.get_interpolated_maxwell_fields(
							e_it.second, E_vals, optical_length_vals);
					o.get_interpolated_maxwell_fields(
							o_it.second, E_vals, optical_length_vals);
					s.get_interpolated_maxwell_fields(
							o_it.first, E_vals, optical_length_vals);

					for(int wave_idx=0; wave_idx<this->N_wavelengths; wave_idx++) {
						for(int pol_idx=0; pol_idx<2; pol_idx++) {
							E_tot = std::complex<double>(0,0);
							for(int ray=0; ray<E_vals[0].size(); ray++)
								E_tot += E_vals[pol_idx][ray]*std::exp(
									I*2.*PI*(optical_length_vals[ray]
										-this->mat_properties.n_host*this->z_step*k) /
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
				MultiDimIndex<2> mesh_idx(
					this->full_horiz_mesh->n_points_per_dim, 0, false);
				mesh_idx = hpx_idx;
				Vector<2,double> screen_pos = 
					this->full_horiz_mesh->origin +
					mesh_idx.get() * this->full_horiz_mesh->cell_lengths;

				int px_idx = hpx_idx + k*this->N_horiz_pixels;
				e.save_interpolated_ray_values(px_idx, screen_pos, this->bulk_extra_data);
				o.save_interpolated_ray_values(px_idx, screen_pos, this->bulk_ordi_data);
			}
		}

		// We evolve the rays to the next z-slice
		std::cout <<
			"	Propagating rays to the next slice" << std::endl << std::endl;

		if(k<this->N_lc_steps-1) {
			#pragma omp parallel for \
				firstprivate(e_seq,o_seq,stepper,trans) private(trans_odes)
			for(int ray_idx=0; ray_idx<N_droplet_rays; ray_idx++) {
				auto& extra_ray_bundle = e.ray_bundles->at(ray_idx);
				auto& ordi_ray_bundle = o.ray_bundles->at(ray_idx);
				int& e_ode_idx = cur_e_ode_idx[ray_idx];
				int& o_ode_idx = cur_o_ode_idx[ray_idx];

				try {
					double z = - this->droplet_radius + (k+1)*this->z_step;

					auto e_full_step = stepper.euler_z_step(
						e_seq.odes[e_ode_idx], extra_ray_bundle,
						s_steps_extra[ray_idx], z);
					auto o_full_step =	stepper.euler_z_step(
						o_seq.odes[o_ode_idx], ordi_ray_bundle,
						s_steps_ordi[ray_idx], z);

					if(!e_full_step && e_ode_idx!=ray_splitting_ode_idx) {
						trans_odes[0] = e_seq.odes[e_ode_idx];
						trans_odes[1] = o_seq.odes[e_ode_idx];
						trans_odes[2] = e_seq.odes[e_ode_idx+1];
						trans_odes[3] = o_seq.odes[e_ode_idx+1];

						trans.update(trans_odes, extra_ray_bundle);
						e_ode_idx++;
					}
					if(!o_full_step && o_ode_idx!=ray_splitting_ode_idx) {
						trans_odes[0] = e_seq.odes[o_ode_idx];
						trans_odes[1] = o_seq.odes[o_ode_idx];
						trans_odes[2] = e_seq.odes[o_ode_idx+1];
						trans_odes[3] = o_seq.odes[o_ode_idx+1];

						trans.update(trans_odes, ordi_ray_bundle);
						o_ode_idx++;
					}
					if(!e_full_step && e_ode_idx==ray_splitting_ode_idx) {
						trans_odes[0] = e_seq.odes[e_ode_idx];
						trans_odes[1] = o_seq.odes[o_ode_idx];
						trans_odes[2] = e_seq.odes[e_ode_idx+1];
						trans_odes[3] = o_seq.odes[o_ode_idx+1];

						trans.update(trans_odes, extra_ray_bundle, ordi_ray_bundle);
						e_ode_idx++;
						o_ode_idx++;
					}
					if(!e_full_step && extra_ray_bundle.cur_pos()(2)<z)
						stepper.euler_z_step(
							e_seq.odes[e_ode_idx], extra_ray_bundle,
							s_steps_extra[ray_idx], z);
					if(!o_full_step && ordi_ray_bundle.cur_pos()(2)<z)
						stepper.euler_z_step(
							o_seq.odes[o_ode_idx], ordi_ray_bundle,
							s_steps_ordi[ray_idx], z);
				}
				catch(...) {
					exc.capture_exception();
				}
			}
			exc.rethrow();

			#pragma omp parallel for firstprivate(s_seq, stepper)
			for(int ray_idx=0; ray_idx<N_iso_rays; ray_idx++) {
				try {
					auto& iso_ray_bundle = s.ray_bundles->at(ray_idx);
					double z = - this->droplet_radius + (k+1)*this->z_step;

					stepper.euler_z_step(
						s_seq.odes[ray_splitting_ode_idx], iso_ray_bundle,
						s_steps_iso[ray_idx], z);
				}
				catch(...) {
					exc.capture_exception();
				}
			}
			exc.rethrow();
		}
		else {
			// We go the upper sample plate
			#pragma omp parallel for \
				firstprivate(e_seq,o_seq,stepper,trans) private(trans_odes)
			for(int ray_idx=0; ray_idx<N_droplet_rays; ray_idx++) {
				auto& extra_ray_bundle = e.ray_bundles->at(ray_idx);
				auto& ordi_ray_bundle = o.ray_bundles->at(ray_idx);
				int& e_ode_idx = cur_e_ode_idx[ray_idx];
				int& o_ode_idx = cur_o_ode_idx[ray_idx];

				while(e_ode_idx!=ray_splitting_ode_idx+3) {
					stepper.euler_step_to_boundary(
						e_seq.odes[e_ode_idx], extra_ray_bundle);

					trans_odes[0] = e_seq.odes[e_ode_idx];
					trans_odes[1] = o_seq.odes[e_ode_idx];
					trans_odes[2] = e_seq.odes[e_ode_idx+1];
					trans_odes[3] = o_seq.odes[e_ode_idx+1];

					trans.update(trans_odes, extra_ray_bundle);
					e_ode_idx++;
				}
				while(o_ode_idx!=ray_splitting_ode_idx+3) {
					stepper.euler_step_to_boundary(
						o_seq.odes[o_ode_idx], ordi_ray_bundle);

					trans_odes[0] = e_seq.odes[o_ode_idx];
					trans_odes[1] = o_seq.odes[o_ode_idx];
					trans_odes[2] = e_seq.odes[o_ode_idx+1];
					trans_odes[3] = o_seq.odes[o_ode_idx+1];

					trans.update(trans_odes, ordi_ray_bundle);
					o_ode_idx++;
				}
			}
			#pragma omp parallel for \
				firstprivate(s_seq,stepper,trans) private(trans_odes)
			for(int ray_idx=0; ray_idx<N_iso_rays; ray_idx++) {
				auto& iso_ray_bundle = s.ray_bundles->at(ray_idx);

				stepper.euler_step_to_boundary(
					s_seq.odes[ray_splitting_ode_idx], iso_ray_bundle);

				trans_odes[0] = s_seq.odes[ray_splitting_ode_idx];
				trans_odes[1] = s_seq.odes[ray_splitting_ode_idx];
				trans_odes[2] = s_seq.odes[ray_splitting_ode_idx+1];
				trans_odes[3] = s_seq.odes[ray_splitting_ode_idx+1];

				trans.update(trans_odes, iso_ray_bundle);
			}
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
		for(int i=0; i<N_droplet_rays; i++) {
			stepper.virtual_z_step(e.ray_bundles->at(i), this->z_foc_1);
			stepper.virtual_z_step(o.ray_bundles->at(i), this->z_foc_1);
		}
		#pragma omp parallel for firstprivate(stepper)
		for(int i=0; i<N_iso_rays; i++)
			stepper.virtual_z_step(s.ray_bundles->at(i), this->z_foc_1);

		std::cout <<
			"	Updating map data..." << std::endl;
		e.update_map_data();
		o.update_map_data();
		s.update_map_data();

		std::cout <<
			"	Inversing screen mapping..." << std::endl;

		std::cout <<
			"		Coarse search of inverses..." << std::endl;

		InverseScreenMap<2> e_inv_map(
			e.target_map, *this->coarse_horiz_mesh, screen_domain,
			this->typical_length, this->hc_parameters);
		InverseScreenMap<2> o_inv_map(
			o.target_map, *this->coarse_horiz_mesh, screen_domain,
			this->typical_length, this->hc_parameters);
		for(int r=0; r<this->N_refinement_cycles; r++) {
			std::cout <<
				"		Refinement of inverses..." << std::endl;
			e_inv_map.refine_data();
			o_inv_map.refine_data();
		}
		auto e_inv_map_data = e_inv_map.get_data();
		auto o_inv_map_data = o_inv_map.get_data();

		// We reconstruct the field values and save them
		std::cout <<
			"	Reconstructing optical fields..." << std::endl;

		#pragma omp parallel for \
			firstprivate(E_tot, e, o, s)
		for(int hpx_idx=0; hpx_idx<this->N_horiz_pixels; hpx_idx++) {
			auto e_it = e_inv_map_data->at(hpx_idx);
			auto o_it = o_inv_map_data->at(hpx_idx);

			std::vector<Vector<3,std::complex<double> > > E_vals[2];
			std::vector<double> optical_length_vals;

			try {
				e.get_interpolated_maxwell_fields(
						e_it.second, E_vals, optical_length_vals);
				o.get_interpolated_maxwell_fields(
						o_it.second, E_vals, optical_length_vals);
				s.get_interpolated_maxwell_fields(
						o_it.first, E_vals, optical_length_vals);

				for(int wave_idx=0; wave_idx<this->N_wavelengths; wave_idx++) {
					for(int pol_idx=0; pol_idx<2; pol_idx++) {
						E_tot = std::complex<double>(0,0);
						for(int ray=0; ray<E_vals[0].size(); ray++)
							E_tot += E_vals[pol_idx][ray]*std::exp(
								I*2.*PI*optical_length_vals[ray] /
								this->wavelengths[wave_idx]);

						for(int comp=0; comp<3; comp++) {
							this->screen_fields_data->E_real[pol_idx][wave_idx]
								->SetComponent(hpx_idx, comp, std::real(E_tot(comp)));
							this->screen_fields_data->E_imag[pol_idx][wave_idx]
								->SetComponent(hpx_idx, comp, std::imag(E_tot(comp)));
						}
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
		for(int i=0; i<N_droplet_rays; i++) {
			stepper.virtual_z_step(e.ray_bundles->at(i), this->lc_thickness/2+1e-4);
			stepper.virtual_z_step(o.ray_bundles->at(i), this->lc_thickness/2+1e-4);
		}
		#pragma omp parallel for firstprivate(stepper)
		for(int i=0; i<N_iso_rays; i++) 
			stepper.virtual_z_step(s.ray_bundles->at(i), this->lc_thickness/2+1e-4);
	}
	

	// If we need ray data on the screen, we propagate rays in the usual
	// way throught the iso layers.
	if(this->screen_ray_output) {
		std::cout <<
			"Interpolating ray data on the focalisation plane..." << std::endl;
		#pragma omp parallel for \
			firstprivate(e_seq, o_seq, stepper, trans) private(trans_odes)
		for(int i=0; i<N_droplet_rays; i++) {
			try {
				auto& extra_ray_bundle = e.ray_bundles->at(i);
				auto& ordi_ray_bundle = o.ray_bundles->at(i);
	
				for(int ode_idx=ray_splitting_ode_idx+3; ode_idx<e_seq.odes.size()-1;
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

				stepper.virtual_z_step(extra_ray_bundle, this->z_foc_2);
				stepper.virtual_z_step(ordi_ray_bundle, this->z_foc_2);
			}
			catch(...) {
				exc.capture_exception();
			}
		}
		exc.rethrow();

		std::cout <<
			"	Updating map data..." << std::endl;
		e.update_map_data();
		o.update_map_data();

		#pragma omp parallel for firstprivate(E_tot, e, o, s)
		for(int hpx_idx=0; hpx_idx<this->N_horiz_pixels; hpx_idx++) {
			MultiDimIndex<2> mesh_idx(
				this->full_horiz_mesh->n_points_per_dim, 0, false);
			mesh_idx = hpx_idx;
			Vector<2,double> screen_pos = 
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
