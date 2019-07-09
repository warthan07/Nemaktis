#include "FresnelTransmission.h"

template <int dim>
void FresnelTransmission<dim>::update(
		std::shared_ptr<ODEFunction<dim> > (&ode_functions)[4],
		RayBundle<dim> &ray_bundle) {

	if(!ode_functions[2]->is_isotropic() && 
			!ode_functions[3]->is_isotropic())
		throw(std::string(
			"Wrong type of ODEs for this fresnel update: the output "
			"medium must correspond to an isotropic phase"));

	// We compute the transmission factors and new polarisations
	for(int i=0; i<2; i++) {
		compute_fresnel_coefs(
			ode_functions, *ray_bundle.rays[0], ray_bundle.cur_pols[i]);

		// If the transmitted rays are evanescent, we set
		// the associated flag, which will allow to skip later
		// propagation steps. 
		if(is_evanescent[2]) {
			ray_bundle.evanescent = true;
			return;
		}

		std::complex<double> a_t1 = fresnel_coefs[2];
		std::complex<double> a_t2 = fresnel_coefs[3];
		double norm = std::sqrt(std::norm(a_t1)+std::norm(a_t2));

		ray_bundle.cur_transmission_factors[i] *=
			norm / std::sqrt(std::abs(ray_bundle.cur_geom_spreading)) /
			ray_bundle.cur_optical_index;
		ray_bundle.cur_pols[i] = (a_t1*e_pol[2] + a_t2*e_pol[3]) / norm;
	}

	// We update the moments and realign the rays to the same optical
	// length
	double max_optical_length = 0;
	for(int i=0; i<dim+1; i++) {
		ode_functions[2]->update_moment(*(ray_bundle.rays[i]));
		max_optical_length =
			std::max(max_optical_length, ray_bundle.rays[i]->cur_optical_length);
	}

	// We add a small quantity to the optical length to void getting
	// stuck on the interface because of rounding error
	max_optical_length += 5*SURFACE_EPS;

	Vector<3,double> pos_derivative, moment_derivative;
	for(int i=0; i<dim+1; i++) {
		if(!ode_functions[2]->def_domain->is_on_surface_embedded(
				ray_bundle.rays[i]->cur_pos))
			throw(std::string("Not on surface"));

		ode_functions[2]->compute_ode_function(
			*(ray_bundle.rays[i]), pos_derivative, moment_derivative);
		ray_bundle.rays[i]->cur_pos +=
			( max_optical_length - ray_bundle.rays[i]->cur_optical_length )
			* pos_derivative;
		ray_bundle.rays[i]->cur_moment +=
			( max_optical_length - ray_bundle.rays[i]->cur_optical_length )
			* moment_derivative;
		ray_bundle.rays[i]->cur_optical_length +=
			( max_optical_length - ray_bundle.rays[i]->cur_optical_length );
	}

	// We update the other ray parameters
	ode_functions[2]->update_optical_index(
		*(ray_bundle.rays[0]), ray_bundle.cur_optical_index);
	ray_bundle.simple_geometrical_spreading_update();
	for(int i=0; i<2; i++)
		ray_bundle.cur_transmission_factors[i] *=
			std::sqrt(std::abs(ray_bundle.cur_geom_spreading)) *
			ray_bundle.cur_optical_index;
}

template <int dim>
void FresnelTransmission<dim>::update(
		std::shared_ptr<ODEFunction<dim> > (&ode_functions)[4],
		RayBundle<dim> &ray_bundle_1,
		RayBundle<dim> &ray_bundle_2) {

	if(ode_functions[2]->is_isotropic() ||
			ode_functions[3]->is_isotropic())
		throw(std::string(
			"Wrong type of ODEs for this fresnel update: the output "
			"medium must correspond to a birefringent phase"));

	// We compute the transmission factor and new polarisation
	for(int i=0; i<2; i++) {
		compute_fresnel_coefs(
			ode_functions, *ray_bundle_1.rays[0], ray_bundle_1.cur_pols[i]);

		if(is_evanescent[2])
			ray_bundle_1.evanescent = true;
		else {
			ray_bundle_1.cur_transmission_factors[i] *=
				2. * fresnel_coefs[2] /
				std::sqrt(std::abs(ray_bundle_1.cur_geom_spreading)) /
				ray_bundle_1.cur_optical_index;
			ray_bundle_1.cur_pols[i] = e_pol[2];
		}

		if(is_evanescent[2])
			ray_bundle_2.evanescent = true;
		else {
			ray_bundle_2.cur_transmission_factors[i] *=
				2. * fresnel_coefs[3] /
				std::sqrt(std::abs(ray_bundle_2.cur_geom_spreading)) /
				ray_bundle_2.cur_optical_index;
			ray_bundle_2.cur_pols[i] = e_pol[3];
		}
	}

	// We update the moments and realign the rays to the same optical
	// length
	Vector<3,double> pos_derivative, moment_derivative;
	if(!is_evanescent[2]) {
		double max_optical_length_1 = 0;

		for(int i=0; i<dim+1; i++) {
			ode_functions[2]->update_moment(*(ray_bundle_1.rays[i]));
			max_optical_length_1 = std::max(
				max_optical_length_1, ray_bundle_1.rays[i]->cur_optical_length);
		}

		// We add a small quantity to the optical length to void getting
		// stuck on the interface because of rounding error
		max_optical_length_1 += 10*SURFACE_EPS;

		for(int i=0; i<dim+1; i++) {
			ode_functions[2]->compute_ode_function(
				*(ray_bundle_1.rays[i]), pos_derivative, moment_derivative);
			ray_bundle_1.rays[i]->cur_pos +=
				(max_optical_length_1 - ray_bundle_1.rays[i]->cur_optical_length)
				* pos_derivative;
			ray_bundle_1.rays[i]->cur_moment +=
				(max_optical_length_1 - ray_bundle_1.rays[i]->cur_optical_length)
				* moment_derivative;
			ray_bundle_1.rays[i]->cur_optical_length +=
				(max_optical_length_1 - ray_bundle_1.rays[i]->cur_optical_length);
		}

		// We update the other ray parameters
		ode_functions[2]->update_optical_index(
				*(ray_bundle_1.rays[0]),
				ray_bundle_1.cur_optical_index);
		ray_bundle_1.simple_geometrical_spreading_update();
		for(int i=0; i<2; i++) 
			ray_bundle_1.cur_transmission_factors[i] *=
				std::sqrt(std::abs(ray_bundle_1.cur_geom_spreading)) *
				ray_bundle_1.cur_optical_index;
	}
	if(!is_evanescent[3]) {
		double max_optical_length_2 = 0;
		for(int i=0; i<dim+1; i++) {
			ode_functions[3]->update_moment(*(ray_bundle_2.rays[i]));
			max_optical_length_2 = std::max(
				max_optical_length_2, ray_bundle_2.rays[i]->cur_optical_length);
		}

		// We add a small quantity to the optical length to void getting
		// stuck on the interface because of rounding error
		max_optical_length_2 += 10*SURFACE_EPS;

		for(int i=0; i<dim+1; i++) {
			ode_functions[3]->compute_ode_function(
				*(ray_bundle_2.rays[i]), pos_derivative, moment_derivative);
			ray_bundle_2.rays[i]->cur_pos +=
				(max_optical_length_2 - ray_bundle_2.rays[i]->cur_optical_length)
				* pos_derivative;
			ray_bundle_2.rays[i]->cur_moment +=
				(max_optical_length_2 - ray_bundle_2.rays[i]->cur_optical_length)
				* moment_derivative;
			ray_bundle_2.rays[i]->cur_optical_length +=
				(max_optical_length_2 - ray_bundle_2.rays[i]->cur_optical_length);
		}

		// We update the other ray parameters
		ode_functions[3]->update_optical_index(
			*(ray_bundle_2.rays[0]),
			ray_bundle_2.cur_optical_index);
		ray_bundle_2.simple_geometrical_spreading_update();
		for(int i=0; i<2; i++) 
			ray_bundle_2.cur_transmission_factors[i] *=
				std::sqrt(std::abs(ray_bundle_2.cur_geom_spreading)) *
				ray_bundle_2.cur_optical_index;
	}
}

template <int dim>
void FresnelTransmission<dim>::compute_fresnel_coefs(
		std::shared_ptr<ODEFunction<dim> > (&ode_functions)[4],
		Ray &ray_i, Vector<3,std::complex<double> > &incident_pol) {

	int i;
	Ray rays[4] = {ray_i, ray_i, ray_i, ray_i};

	Vector<3,double> normal;
	if(!ode_functions[0]->def_domain->get_normal_embedded(
			ray_i.cur_pos, normal))
		throw(std::string(
			"Error: you are trying to apply the fresnel boundary conditions "
			"but the ray bundle is not on a surface"));

	e_i = incident_pol;
	b_i = ray_i.cur_moment ^ e_i;

	t_s = ray_i.cur_moment ^ normal;
	if(t_s.norm()<1e-10) {
		t_s = {0., normal(2), -normal(1)};
		if(t_s.norm()<1e-10)
			t_s = {-normal(1), normal(0), 0.};
	}
	t_s /= t_s.norm();
	t_p = normal ^ t_s;
	t_p /= t_p.norm();

	// We initialize the polarisation vectors with s and
	// p polarizations if the medium is isotropic(), or with 
	// eigenvectors polarisation if the medium is anisotropic().
	Vector<3,std::complex<double> > moment;
	for(i=0; i<2; i++) {
		ode_functions[2*i]->update_moment(rays[2*i]);
		is_evanescent[2*i] = rays[2*i].is_evanescent;

		ode_functions[2*i+1]->update_moment(rays[2*i+1]);
		is_evanescent[2*i+1] = rays[2*i+1].is_evanescent;

		if(ode_functions[2*i]->is_isotropic() &&
				ode_functions[2*i+1]->is_isotropic()) {

			if(is_evanescent[2*i])
				moment = rays[2*i].complex_moment;
			else
				moment = rays[2*i].cur_moment;

			e_pol[2*i] = moment ^ normal;
			if(e_pol[2*i].norm()<1e-10) {
				e_pol[2*i] = {0., normal(2), -normal(1)};
				if(e_pol[2*i].norm()<1e-10)
					e_pol[2*i] = {-normal(1), normal(0), 0.};
			}
			e_pol[2*i] /= e_pol[2*i].norm();
			b_pol[2*i] = moment ^ e_pol[2*i];

			if(is_evanescent[2*i+1])
				moment = rays[2*i+1].complex_moment;
			else
				moment = rays[2*i+1].cur_moment;

			e_pol[2*i+1] = e_pol[2*i] ^ moment;
			e_pol[2*i+1] /= e_pol[2*i+1].norm();
			b_pol[2*i+1] = moment ^ e_pol[2*i+1];
		}
		else if(!ode_functions[2*i]->is_isotropic() && 
				!ode_functions[2*i+1]->is_isotropic()) {

			e_pol[2*i] = e_i;
			ode_functions[2*i]->compute_polarisation(rays[2*i], e_pol[2*i], true);
			b_pol[2*i] = rays[2*i].cur_moment ^ e_pol[2*i];		


			e_pol[2*i+1] = b_i/b_i.norm();
			ode_functions[2*i+1]->compute_polarisation(rays[2*i+1], e_pol[2*i+1], true);
			b_pol[2*i+1] = rays[2*i+1].cur_moment ^ e_pol[2*i+1];		
		}
		else
			throw(std::string(
				"Wrong ODEFunction array: the first (last) two components\n"
				"must have the same is_isotropic() value (FresnelTransmission)"));
	}

	// We assemble the Fresnel system
	fresnel_matrix <<
		(t_s,e_pol[0]), (t_s,e_pol[1]), (t_s,e_pol[2]), (t_s,e_pol[3]),
		(t_p,e_pol[0]), (t_p,e_pol[1]), (t_p,e_pol[2]), (t_p,e_pol[3]),
		(t_s,b_pol[0]), (t_s,b_pol[1]), (t_s,b_pol[2]), (t_s,b_pol[3]),
		(t_p,b_pol[0]), (t_p,b_pol[1]), (t_p,b_pol[2]), (t_p,b_pol[3]);
	fresnel_rhs <<
		(t_s,e_i),
		(t_p,e_i),
		(t_s,b_i),
		(t_p,b_i);

	// We solve the Fresnel system using the Eigen HouseholderQR decomposition
	// algorithm
	fresnel_coefs = fresnel_matrix.householderQr().solve(fresnel_rhs);
}


template class FresnelTransmission<1>;
template class FresnelTransmission<2>;
template class FresnelTransmission<3>;
