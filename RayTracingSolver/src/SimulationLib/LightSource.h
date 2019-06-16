#ifndef LIGHTSOURCE_H
#define LIGHTSOURCE_H

#include <memory>

#include "CartesianMesh.h"
#include "DefinitionDomain.h"
#include "RayBundle.h"

template <int dim>
class LightSource {
public:
	LightSource(
		json settings,
		const std::shared_ptr<CartesianMesh<dim> > &source_mesh,
		const std::shared_ptr<DefinitionDomain<dim> > &source_domain);

	const std::shared_ptr<CartesianMesh<dim> >& get_source_mesh() const;
	const std::shared_ptr<DefinitionDomain<dim> >& get_source_domain() const;

	void init_rays(
		std::shared_ptr<std::vector<RayBundle<dim+1> > > &ray_bundles,
		std::vector<int> &ray_to_mesh_indices,
		double z_init, double init_ampl) const;

private:
	std::shared_ptr<CartesianMesh<dim> > source_mesh;
	std::shared_ptr<DefinitionDomain<dim> > source_domain;

	Vector<3,double> pol_basis[2];
	double eps_r;
	double mu;
};

#endif
