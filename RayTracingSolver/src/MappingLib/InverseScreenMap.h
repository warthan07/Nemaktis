#ifndef INVERSESCREENMAP_H
#define INVERSESCREENMAP_H

#include "CubicInterpolatedMapping.h"
#include "HCInverseSolver.h"

template <int dim>
class InverseScreenMap {
public:
	InverseScreenMap(
		CubicInterpolatedMapping<dim,dim,double> &screen_map,
		const CartesianMesh<dim> &coarse_mesh,
		std::shared_ptr<DefinitionDomain<dim> > &screen_domain,
		double typical_length,
		json settings);

	const std::vector<DataPair<dim> >* get_data() const {
		return &coarse_data;
	}
	void refine_data();

	const Vector<dim,unsigned long>& n_points_per_dim() const {
		return coarse_mesh.n_points_per_dim;
	}

private:
	CubicInterpolatedMapping<dim,dim,double> screen_map;

	CartesianMesh<dim> coarse_mesh, fine_mesh;
	std::vector<DataPair<dim> > coarse_data, fine_data;
	std::shared_ptr<DefinitionDomain<dim> > screen_domain;

	double typical_length;
	json coarse_hc_parameters, refinement_hc_parameters;
};

#endif
