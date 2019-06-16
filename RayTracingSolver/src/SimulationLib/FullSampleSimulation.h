#ifndef FULLSAMPLESIMULATION_H
#define FULLSAMPLESIMULATION_H

#include "Simulation.h"
#include "RayStepper.h"
#include "FresnelTransmission.h"

template <int dim>
class FullSampleSimulation : public Simulation<dim> {
public:
	FullSampleSimulation(
		json settings, double lc_thickness, unsigned int N_lc_steps,
		std::shared_ptr<CubicInterpolatedMapping<dim,3,double> > &n_field);

	virtual void run();

private:
	std::shared_ptr<LightSource<dim-1> > source;

	std::shared_ptr<RayFamily<dim> > extra_ray_family, ordi_ray_family;
	ODESequence<dim> extra_ray_ode_seq, ordi_ray_ode_seq;
	std::shared_ptr<ODEFunction<dim> > extra_ode, ordi_ode;

	int N_rays;
	int ray_splitting_ode_idx;
};

#endif
