#include "RayStepper.h"

template <int dim>
bool RayStepper<dim>::euler_step(
		std::shared_ptr<ODEFunction<dim> > &ode_function,
		RayBundle<dim> &ray_bundle,
		double &step_size) {
	
	return adaptative_euler_step(
		ode_function, ray_bundle, step_size, false);
}

template <int dim>
bool RayStepper<dim>::euler_z_step(
		std::shared_ptr<ODEFunction<dim> > &ode_function,
		RayBundle<dim> &ray_bundle,
		double &step_size,
		double target_z) {

	// We apply euler steps until the z_target plane is crossed
	// by the first ray of the bundle or until a boundary is crossed
	bool reached_boundary = false;
	while(ray_bundle.rays[0]->cur_pos(2)<target_z) {
		if(!adaptative_euler_step(ode_function, ray_bundle, step_size, true)) {
			reached_boundary = true;
			break;
		}
	}

	// If a boundary_step was deemed necessary in the previous loop,
	// we check if the target plane will be crossed or not during this
	// boundary step. If no, we can apply the full boundary step and
	// return false to indicate that this z_step is incomplete. If yes,
	// we can go instead to the z_target plane and return true to
	// indicate that the z_step is complete.
	if(reached_boundary) {
		Vector<3,double> new_pos;
		ode_function->def_domain->compute_intersection_embedded(
			ray_bundle.rays[0]->cur_pos, pos_derivative[0], new_pos);

		// The shift with fd_step() is necessary in case target_z is
		// aligned with the interface (it is simpler to just shift the
		// target plane and apply boundary step and fresnel condition
		// as usual).
		double shift = 200*ray_bundle.get_fd_step();
		if(new_pos(2)<target_z-shift) {
			boundary_step(ode_function, ray_bundle);
			return false;
		}	
		else if(new_pos(2)<target_z+shift) {
			bulk_step(
				ode_function, ray_bundle,
				(target_z-shift-ray_bundle.rays[0]->cur_pos(2))
					/ pos_derivative[0](2));
			return true;
		}
		else {
			bulk_step(
				ode_function, ray_bundle,
				(target_z-ray_bundle.rays[0]->cur_pos(2)) / pos_derivative[0](2));
			return true;
		}
	}
	else {
		// We go to the z_target plane and return true to indicate that
		// the z_step is complete.
		bulk_step(
			ode_function, ray_bundle,
			(target_z - ray_bundle.rays[0]->cur_pos(2)) / pos_derivative[0](2));
		return true;
	}
}

template <int dim>
void RayStepper<dim>::euler_step_to_boundary(
		std::shared_ptr<ODEFunction<dim> > &ode_function,
		RayBundle<dim> &ray_bundle) {

	for(int i=0; i<dim+1; i++) {
		Ray &ray = (*ray_bundle.rays[i]);
		ode_function->compute_ode_function(
			ray, pos_derivative[i], moment_derivative[i]);
	}
	boundary_step(ode_function, ray_bundle);
}

template <int dim>
void RayStepper<dim>::virtual_z_step(
		RayBundle<dim> &ray_bundle, double target_z) {

	double eps_r = std::pow(ray_bundle.cur_optical_index,2.);
	double optical_length_step =
		( target_z - ray_bundle.rays[0]->cur_pos(2) ) *
		eps_r / ray_bundle.rays[0]->cur_moment(2);

	Vector<3,double> pos_derivatives[dim+1];
	for(int i=0; i<dim+1; i++) {
		Ray &ray = *(ray_bundle.rays[i]);

		pos_derivatives[i] = ray.cur_moment /eps_r;
		ray.cur_pos += optical_length_step * pos_derivatives[i];
		ray.cur_optical_length += optical_length_step;
	}

	ray_bundle.exact_geometrical_spreading_update(
		pos_derivatives, optical_length_step);
}

template <int dim>
bool RayStepper<dim>::adaptative_euler_step(
		std::shared_ptr<ODEFunction<dim> > &ode_function,
		RayBundle<dim> &ray_bundle,
		double &step_size,
		bool do_not_update_if_boundary_step) {

	// We compute the ode functions for the position and moment
	// and check if the specified step length will yield a new position
	// outside the current definition domain. 
	Vector<3,double> new_pos;
	bool reached_boundary = false;
	for(int i=0; i<dim+1; i++) {
		Ray &ray = (*ray_bundle.rays[i]);
		ode_function->compute_ode_function(
			ray, pos_derivative[i], moment_derivative[i]);

		new_pos = ray.cur_pos + step_size * pos_derivative[i];
		if(!ode_function->def_domain->is_inside_embedded(
					new_pos,-SURFACE_EPS))

			reached_boundary = true;
	}

	if(!reached_boundary) {
		// Adaptive change of the step size based on the new value of
		// the hamiltonian after the euler step. We store the boolean
		// returned by update_step_size to determine if we need to
		// increase the step size after updating the state of the ray
		// bundle.
		bool increase_step = update_step_size(
			ode_function, ray_bundle, step_size);

		// We update the ray bundle based on the calculated step size
		bulk_step(ode_function, ray_bundle, step_size);

		// We increase the step if the error was lower than half the
		// threshold
		if(increase_step)
			step_size *= 2;

		return true;
	}
	else {
		// Adaptive change of the step size based on the new value of
		// the hamiltonian after the euler step. If the step size can be
		// increased (first condition), this means that we can safely
		// use a boundary step. If the step size was decreased (second
		// condition), we still need a few bulk step before going to the
		// boundary.
		if(update_step_size(ode_function, ray_bundle, step_size)) {
			if(!do_not_update_if_boundary_step)
				boundary_step(ode_function, ray_bundle);
			return false;
		}
		else {
			bulk_step(ode_function, ray_bundle, step_size);
			return true;
		}
	}
}

template <int dim>
bool RayStepper<dim>::update_step_size(
		std::shared_ptr<ODEFunction<dim> > &ode_function,
		RayBundle<dim> &ray_bundle,
		double &step_size) {

	// We decrease the step size until the new positions of the rays are
	// inside the domain
	Vector<3,double> pos;
	for(int i=0; i<dim+1; i++) {
		pos = ray_bundle.rays[i]->cur_pos + step_size * pos_derivative[i];
		while(!ode_function->def_domain->is_inside_embedded(pos,-SURFACE_EPS)) {
			step_size /= 2;
			pos -= step_size * pos_derivative[i];
		}
	}

	// We decrease the step size until the new state is associated with
	// an error below the threshold for the 0-th ray of the bundle.
	Ray ray = *(ray_bundle.rays[0]);
	ray.cur_pos += step_size * pos_derivative[0];
	ray.cur_moment += step_size * moment_derivative[0];

	bool increase_step = false;
	double error = std::abs(2.*ode_function->hamiltonian_value(ray)-1.);
	if(error>threshold)
		while(error>threshold) {
			step_size /= 2;
			ray.cur_pos -= step_size * pos_derivative[0];
			ray.cur_moment -= step_size * moment_derivative[0];
			error = std::abs(2.*ode_function->hamiltonian_value(ray)-1.);
		}
	else if(error<threshold/2 && 2*step_size <= max_step_size)
		return true; // indicate that we can increase the step size 

	return false; // indicate that we do not need to change the step size
}

template <int dim>
void RayStepper<dim>::bulk_step(
		std::shared_ptr<ODEFunction<dim> > &ode_function,
		RayBundle<dim> &ray_bundle,
		double step_size) {

	// We update the rays in the bundle
	for(int i=0; i<dim+1; i++) {
		Ray &ray = *(ray_bundle.rays[i]);
		ray.cur_pos += step_size * pos_derivative[i];
	 	ray.cur_moment += step_size * moment_derivative[i];
		ray.cur_optical_length += step_size;

		// We apply a simple correction step to preserve the energy
		ray.cur_moment *=
			1. / std::sqrt(2.*ode_function->hamiltonian_value(ray));
	}

	// We update the geometrical spreading
	ray_bundle.exact_geometrical_spreading_update(
		pos_derivative, step_size);

	// Finally, we update the polarisation
	for(int i=0; i<2; i++)
		ode_function->compute_polarisation(
			*(ray_bundle.rays[0]), ray_bundle.cur_pols[i]);
}

template <int dim>
void RayStepper<dim>::boundary_step(
		std::shared_ptr<ODEFunction<dim> > &ode_function,
		RayBundle<dim> &ray_bundle) {

	// We compute the minimum step size allowing at least one ray to
	// reach the boundary.
	int i;
	double min_step_size;

	Vector<3,double> new_pos[dim+1];
	Vector<3,double> new_moment[dim+1];

	step_sizes[0] =
		ode_function->def_domain->compute_intersection_embedded(
			ray_bundle.rays[0]->cur_pos, pos_derivative[0], new_pos[0]);
	if(step_sizes[0]<0) {
		throw(std::string(
			"Could not find an intersection between the ray bundle and "
			"the specified domain. Check your definition of the domains."));
	}

	min_step_size = step_sizes[0];
	for(i=1; i<dim+1; i++) {
		step_sizes[i] =
			ode_function->def_domain->compute_intersection_embedded(
				ray_bundle.rays[i]->cur_pos, pos_derivative[i], new_pos[i]);
		if(step_sizes[i]<0) {
			throw(std::string(
				"Could not find an intersection between the ray bundle and "
				"the specified domain. Check your definition of the domains."));
		}
		min_step_size = std::min(min_step_size, step_sizes[i]);
	}

	// We evolve the rays using the minimum step size and update
	// the RayBundle parameters
	for(i=0; i<dim+1; i++) {
		Ray &ray = *(ray_bundle.rays[i]);
		new_moment[i] = ray.cur_moment + step_sizes[i] * moment_derivative[i];

		ray.cur_pos += min_step_size * pos_derivative[i];
		ray.cur_moment += min_step_size * moment_derivative[i];
	}
	ray_bundle.exact_geometrical_spreading_update(
		pos_derivative, min_step_size);

	for(int i=0; i<2; i++)
		ode_function->compute_polarisation(
			*(ray_bundle.rays[0]), ray_bundle.cur_pols[i]);
	ode_function->update_optical_index(
		*(ray_bundle.rays[0]),
		ray_bundle.cur_optical_index);

	// We finish evolving the rays to the interface using the full
	// step sizes.
	for(i=0; i<dim+1; i++) {
		Ray &ray = *(ray_bundle.rays[i]);
		ray.cur_pos = new_pos[i];
		ray.cur_moment = new_moment[i];
		ray.cur_optical_length += step_sizes[i];

		// We apply a simple correction step to preserve the energy
		ray.cur_moment *= 1. / std::sqrt(2.*ode_function->hamiltonian_value(ray));
	}
}

template class RayStepper<2>;
template class RayStepper<3>;
