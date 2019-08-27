#include "RayFamily.h"

template <int dim>
RayFamily<dim>::RayFamily(LightSource<dim-1> &source, double z_init, double init_ampl) : 
	target_map_data(std::make_shared<std::vector<Vector<dim-1,double> > >(
		source.get_source_mesh()->n_points)),
	E_maps_data{
		std::make_shared<std::vector<Vector<3,std::complex<double> > > >(
			source.get_source_mesh()->n_points),
		std::make_shared<std::vector<Vector<3,std::complex<double> > > >(
			source.get_source_mesh()->n_points)},
	moment_map_data(std::make_shared<std::vector<Vector<3,double> > >(
		source.get_source_mesh()->n_points)),
	optical_length_map_data(std::make_shared<std::vector<Vector<1,double> > >(
		source.get_source_mesh()->n_points)),
	NA_limited_map_data(std::make_shared<std::vector<Vector<1,double> > >(
		source.get_source_mesh()->n_points)),
	evanescent_map_data(std::make_shared<std::vector<Vector<1,double> > >(
		source.get_source_mesh()->n_points)),
	target_map(
		target_map_data, source.get_source_mesh(), source.get_source_domain()),
	E_maps{
		{E_maps_data[0], source.get_source_mesh(), source.get_source_domain()},
		{E_maps_data[1], source.get_source_mesh(), source.get_source_domain()}},
	moment_map(
		moment_map_data, source.get_source_mesh(), source.get_source_domain()),
	optical_length_map(
		optical_length_map_data, source.get_source_mesh(), source.get_source_domain()),
	NA_limited_map(
		NA_limited_map_data, source.get_source_mesh(), source.get_source_domain()),
	evanescent_map(
		evanescent_map_data, source.get_source_mesh(), source.get_source_domain()) {

	source.init_rays(ray_bundles, ray_to_mesh_indices, z_init, init_ampl);
}

template <int dim>
void RayFamily<dim>::update_map_data() {

	#pragma omp parallel for
	for(int i=0; i<ray_bundles->size(); i++) {
		Vector<3, std::complex<double> > E;
		for(int pol_idx=0; pol_idx<2; pol_idx++) {
			ray_bundles->at(i).compute_E_field(E, pol_idx);
			E_maps_data[pol_idx]->at(ray_to_mesh_indices[i]) = E;
		}

		optical_length_map_data->at(ray_to_mesh_indices[i]) =
			ray_bundles->at(i).optical_length();

		moment_map_data->at(ray_to_mesh_indices[i]) =
			ray_bundles->at(i).cur_moment();

		Vector<3,double> pos = ray_bundles->at(i).cur_pos();
		for(int j=0; j<dim-1; j++)
			target_map_data->at(ray_to_mesh_indices[i])(j) = pos(j+3-dim);

		if(ray_bundles->at(i).is_NA_limited())
			NA_limited_map_data->at(ray_to_mesh_indices[i])(0) = 1;
		else
			NA_limited_map_data->at(ray_to_mesh_indices[i])(0) = 0;

		if(ray_bundles->at(i).is_evanescent())
			evanescent_map_data->at(ray_to_mesh_indices[i])(0) = 1;
		else
			evanescent_map_data->at(ray_to_mesh_indices[i])(0) = 0;
	}
	
	for(int pol_idx=0; pol_idx<2; pol_idx++)
		E_maps[pol_idx].extrapolate_data();
	optical_length_map.extrapolate_data();
	moment_map.extrapolate_data();
	target_map.extrapolate_data();
	NA_limited_map.extrapolate_data();
	evanescent_map.extrapolate_data();
}

template <int dim>
void RayFamily<dim>::get_interpolated_maxwell_fields(
		const Vector<dim-1,double> &pos_val,
		std::vector<Vector<3,std::complex<double> > > (&E_vals)[2],
		std::vector<double> &optical_length_vals) {

	if(NA_limited_map.get_value(pos_val,NA_limited) && 
			std::abs(NA_limited(0))<0.5 &&
		evanescent_map.get_value(pos_val,evanescent) && 
			std::abs(evanescent(0))<0.5) {

		optical_length_map.get_value(pos_val,optical_length);
		optical_length_vals.push_back(optical_length(0));

		for(int pol_idx=0; pol_idx<2; pol_idx++) {
			E_maps[pol_idx].get_value(pos_val,E);
			E_vals[pol_idx].push_back(E);
		}
	}
}

template <int dim>
void RayFamily<dim>::get_interpolated_maxwell_fields(
		const std::vector<Vector<dim-1,double> > &pos_vals,
		std::vector<Vector<3,std::complex<double> > > (&E_vals)[2],
		std::vector<double> &optical_length_vals) {

	for(auto x : pos_vals)
		get_interpolated_maxwell_fields(
			x, E_vals, optical_length_vals);
}

template <int dim>
void RayFamily<dim>::save_interpolated_ray_values(
		unsigned int px_idx, 
		const Vector<dim-1,double> &pos_val,
		std::shared_ptr<RayData> &ray_data) {

	if(evanescent_map.get_value(pos_val,evanescent) && 
			std::abs(evanescent(0))<0.5) {

		// We assume that the imaginary part in E_maps is 0
		for(int pol_idx=0; pol_idx<2; pol_idx++) {
			E_maps[pol_idx].get_value(pos_val,E);
			for(int comp=0; comp<3; comp++)
				ray_data->ampl[pol_idx]->SetComponent(
					px_idx, comp, std::real(E(comp)));
		}

		moment_map.get_value(pos_val,moment);
		for(int comp=0; comp<3; comp++)
			ray_data->moment->SetComponent(px_idx, comp, moment(comp));

		optical_length_map.get_value(pos_val,optical_length);
		ray_data->opt_length->SetComponent(px_idx, 0, optical_length(0));

		target_map.get_value(pos_val, r);
		switch(dim) {
			case 2:
				ray_data->deflection->SetComponent(px_idx, 0, 0);
				ray_data->deflection->SetComponent(px_idx, 1, r(0)-pos_val(0));
				break;
			case 3:
				ray_data->deflection->SetComponent(px_idx, 0, r(0)-pos_val(0));
				ray_data->deflection->SetComponent(px_idx, 1, r(1)-pos_val(1));
				break;
		}
	}
	else {
		for(int pol_idx=0; pol_idx<2; pol_idx++)
			for(int comp=0; comp<3; comp++)
				ray_data->ampl[pol_idx]->SetComponent(px_idx, comp, 0);
		for(int comp=0; comp<2; comp++)
			ray_data->deflection->SetComponent(px_idx, comp, 0);
		ray_data->opt_length->SetComponent(px_idx, 0, 0);
	}
}

template class RayFamily<2>;
template class RayFamily<3>;
