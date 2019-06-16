#ifndef DROPLETSAMPLESIMULATION_H
#define DROPLETSAMPLESIMULATION_H

#include "Simulation.h"
#include "RayStepper.h"
#include "FresnelTransmission.h"

class DropletSampleSimulation : public Simulation<3> {
public:
	DropletSampleSimulation(
		json settings, double lc_thickness, unsigned int N_lc_steps,
		std::shared_ptr<CubicInterpolatedMapping<3,3,double> > &n_field);

	virtual void run();

private:
	std::shared_ptr<LightSource<2> > droplet_source, iso_source;

	std::shared_ptr<RayFamily<3> >
		extra_ray_family, ordi_ray_family, iso_ray_family;
	ODESequence<3> extra_ray_ode_seq, ordi_ray_ode_seq, iso_ray_ode_seq;
	std::shared_ptr<ODEFunction<3> > extra_ode, ordi_ode, iso_ode;

	double droplet_radius;
	int N_droplet_rays, N_iso_rays;
	int ray_splitting_ode_idx;
};

#endif
