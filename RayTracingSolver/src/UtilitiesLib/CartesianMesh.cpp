#include "CartesianMesh.h"

template <int dim>
CartesianMesh<dim>::CartesianMesh(
		const Vector<dim,double> &mesh_origin,
		const Vector<dim,double> &mesh_lengths,
		const Vector<dim,long> &n_points_per_dim,
		bool last_idx_fastest) :
	origin(mesh_origin),
	lengths(mesh_lengths),
	n_points_per_dim(n_points_per_dim),
	last_idx_fastest(last_idx_fastest) {

	int d;

	n_points=1;
	for(d=0; d<dim; ++d) {
		cell_lengths(d) = lengths(d) / (n_points_per_dim(d)-1);
		n_points *= n_points_per_dim(d);
	}

	// Total number of points on a d-dimensional cell
	n_points_per_cell = 1;
	for(d=0; d<dim; d++)
		n_points_per_cell *= 2;

	if(last_idx_fastest) {
		flatten_weights(dim-1) = 1;
		for(d=dim-2; d>=0; --d)
			flatten_weights(d) = flatten_weights(d+1)*n_points_per_dim(d+1);
	}
	else {
		flatten_weights(0) = 1;
		for(d=1; d<dim; ++d)
			flatten_weights(d) = flatten_weights(d-1)*n_points_per_dim(d-1);
	}
}

template <int dim>
void CartesianMesh<dim>::refine() {

	int d;

	cell_lengths /= 2.;
	n_points = 1;
	for(d=0; d<dim; ++d) {
		n_points_per_dim(d) = 2*n_points_per_dim(d)-1;
		n_points *= n_points_per_dim(d);
	}

	if(last_idx_fastest) {
		flatten_weights(dim-1) = 1;
		for(d=dim-2; d>=0; --d)
			flatten_weights(d) = flatten_weights(d+1)*n_points_per_dim(d+1);
	}
	else {
		flatten_weights(0) = 1;
		for(d=1; d<dim; ++d)
			flatten_weights(d) = flatten_weights(d-1)*n_points_per_dim(d-1);
	}
}

template class CartesianMesh<1>;
template class CartesianMesh<2>;
template class CartesianMesh<3>;
