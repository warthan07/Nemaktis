#ifndef RAYSTEPPER_H
#define RAYSTEPPER_H

#include "ODEFunction.h"

template <int dim>
class RayStepper {
public:
	RayStepper(double max_step_size, double threshold = 1e-4) :
		max_step_size(max_step_size),
		threshold(threshold) {}

	bool euler_step(
		std::shared_ptr<ODEFunction<dim> > &ode_function,
		RayBundle<dim> &ray_bundle,
		double &step_size);

	bool euler_z_step(
		std::shared_ptr<ODEFunction<dim> > &ode_function,
		RayBundle<dim> &ray_bundle,
		double &step_size,
		double target_z);

	void euler_step_to_boundary(
		std::shared_ptr<ODEFunction<dim> > &ode_function,
		RayBundle<dim> &ray_bundle);

	void virtual_z_step(
		RayBundle<dim> &ray_bundle,
		double target_z);

protected:
	bool adaptative_euler_step(
		std::shared_ptr<ODEFunction<dim> > &ode_function,
		RayBundle<dim> &ray_bundle,
		double &step_size,
		bool do_not_update_if_boundary_step);
	bool update_step_size(
		std::shared_ptr<ODEFunction<dim> > &ode_function,
		RayBundle<dim> &ray_bundle,
		double &step_size);
	void bulk_step(
		std::shared_ptr<ODEFunction<dim> > &ode_function,
		RayBundle<dim> &ray_bundle,
		double step_size);
	void boundary_step(
		std::shared_ptr<ODEFunction<dim> > &ode_function,
		RayBundle<dim> &ray_bundle);

	Vector<3,double> pos_derivative[dim+1];
	Vector<3,double> moment_derivative[dim+1];

	double step_sizes[dim+1];
	double max_step_size;

	double threshold;
};

#endif
