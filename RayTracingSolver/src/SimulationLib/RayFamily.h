#ifndef RAYFAMILY_H
#define RAYFAMILY_H

#include "LightSource.h"
#include "CubicInterpolatedMapping.h"
#include "VTIData.h"

template <int dim>
class RayFamily {
public:
	RayFamily(LightSource<dim-1> &source, double z_init, double init_ampl);

	void update_map_data();

	void get_interpolated_maxwell_fields(
		const Vector<dim-1,double> &pos_val,
		std::vector<Vector<3,std::complex<double> > > (&E_vals)[2],
		std::vector<double> &optical_length_vals);
	void get_interpolated_maxwell_fields(
		const std::vector<Vector<dim-1,double> > &pos_vals,
		std::vector<Vector<3,std::complex<double> > > (&E_vals)[2],
		std::vector<double> &optical_length_vals);

	void save_interpolated_ray_values(
		unsigned int px_idx, 
		const Vector<dim-1,double> &pos_val,
		std::shared_ptr<RayData> &ray_data);

	std::shared_ptr<std::vector<RayBundle<dim> > > ray_bundles;
	std::vector<int> ray_to_mesh_indices;

	std::shared_ptr<std::vector<Vector<dim-1,double> > > target_map_data;
	std::shared_ptr<std::vector<Vector<3,std::complex<double> > > > E_maps_data[2];
	std::shared_ptr<std::vector<Vector<3,double> > > moment_map_data;
	std::shared_ptr<std::vector<Vector<1,double> > > optical_length_map_data;
	std::shared_ptr<std::vector<Vector<1,double> > > NA_limited_map_data;
	std::shared_ptr<std::vector<Vector<1,double> > > evanescent_map_data;

	CubicInterpolatedMapping<dim-1,dim-1,double> target_map;
	CubicInterpolatedMapping<dim-1,3,std::complex<double> > E_maps[2];
	CubicInterpolatedMapping<dim-1,3,double> moment_map;
	CubicInterpolatedMapping<dim-1,1,double> optical_length_map;
	CubicInterpolatedMapping<dim-1,1,double> NA_limited_map;
	CubicInterpolatedMapping<dim-1,1,double> evanescent_map;

private:
	Vector<3, std::complex<double> > E;
	Vector<3,double> moment;
	Vector<dim-1,double> r;
	Vector<1,double> optical_length, NA_limited, evanescent;
};

#endif
