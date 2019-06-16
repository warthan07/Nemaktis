#include "LightSource.h"

template <int dim>
LightSource<dim>::LightSource(
		json settings,
		const std::shared_ptr<CartesianMesh<dim> > &source_mesh,
		const std::shared_ptr<DefinitionDomain<dim> > &source_domain) :
	source_mesh(source_mesh),
	source_domain(source_domain) {

	pol_basis[0] = {1,0,0};
	pol_basis[1] = {0,1,0};

	eps_r = std::pow(
		settings.at("Material properties").at(
			"External medium refractive index"), 2.);
	mu = source_mesh->lengths.norm() / 1e8;
}

template <int dim>
const std::shared_ptr<CartesianMesh<dim> >&
LightSource<dim>::get_source_mesh() const {
	return source_mesh;
}

template <int dim>
const std::shared_ptr<DefinitionDomain<dim> >&
LightSource<dim>::get_source_domain() const {
	return source_domain;
}

template <int dim>
void LightSource<dim>::init_rays(
		std::shared_ptr<std::vector<RayBundle<dim+1> > > &ray_bundles,
		std::vector<int> &ray_to_mesh_indices,
		double z_init, double init_ampl) const {

	MultiDimIndex<dim> idx(source_mesh->n_points_per_dim);

	Vector<3,double> p({0., 0., std::sqrt(eps_r)});
	Vector<3,double> x({0., 0., z_init});
	Vector<dim,double> xs;

	ray_bundles = std::make_shared<std::vector<RayBundle<dim+1> > >();
	ray_to_mesh_indices.clear();

	for(; idx.valid();  idx++) {
		xs = source_mesh->origin + idx.get()*source_mesh->cell_lengths;
		if(source_domain->is_inside(xs)) {
			for(int i=0; i<dim; i++)
				x(i+2-dim) = xs(i);

			auto ray_bundle = RayBundle<dim+1>(
				x, p, pol_basis, std::sqrt(eps_r), init_ampl,  mu);
			ray_bundles->push_back(ray_bundle);
			ray_to_mesh_indices.push_back(idx());
		}
	}
}

template class LightSource<1>;
template class LightSource<2>;
