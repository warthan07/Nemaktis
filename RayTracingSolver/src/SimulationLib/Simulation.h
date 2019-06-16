#ifndef SIMULATION_H
#define SIMULATION_H

#include "RayFamily.h"
#include "LightSource.h"
#include "ODEFunction.h"
#include "CubicInterpolatedMapping.h"
#include "VTIData.h"
#include "MaterialProperties.h"

template <int dim>
class Simulation {
public:
	Simulation(
		json settings, double lc_thickness, unsigned int N_lc_steps,
		std::shared_ptr<CubicInterpolatedMapping<dim,3,double> > &n_field);

	virtual void run() = 0;

protected:
	void apply_fourier_iso_layer_filter();

	std::shared_ptr<CartesianMesh<dim-1> > coarse_horiz_mesh, full_horiz_mesh;
	int N_refinement_cycles;
	json hc_parameters;

	json viz_settings;
	std::string bulk_basename, screen_basename, basedir;
	bool bulk_ray_output, bulk_fields_output, screen_ray_output, screen_fields_output;
	std::shared_ptr<RayData> bulk_extra_data, screen_extra_data;
	std::shared_ptr<RayData> bulk_ordi_data, screen_ordi_data;
	std::shared_ptr<FieldsData> bulk_fields_data, screen_fields_data;

	double z_step;
	unsigned long N_lc_steps, N_horiz_pixels;

	std::vector<double> wavelengths;
	unsigned int N_wavelengths;
	double numerical_aperture;

	MaterialProperties mat_properties;
	double typical_length;
	double lc_thickness, lower_iso_thickness, upper_iso_thickness;
	double z_foc_1, z_foc_2;
	std::vector<std::shared_ptr<ODEFunction<dim> > >
		lower_iso_odes, upper_iso_odes;

	std::shared_ptr<CubicInterpolatedMapping<dim,3,double> > n_field;
};

#endif
