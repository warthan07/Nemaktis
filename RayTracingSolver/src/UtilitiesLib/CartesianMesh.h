#ifndef CARTESIANMESH_H
#define CARTESIANMESH_H

#include "linalg.h"

/**
 * A container class which allows to fully caracterize a cartesian mesh
 * in dimension \f$d_1=\f$<tt>dim</tt>.
 * 
 * The nodes \f$\vec{x}\left[i_1,...i_{d_1}\right]\f$ of the cartesian
 * mesh are defined as:
 * \f[
 * 		\vec{x}\left[i_1,...i_{d_1}\right] = \vec{x}_0 + \sum_{d=1}^{d_1}
 * 			\frac{i_d-1}{N_d-1}\,L_d\,\vec{e}_d,
 * \f]
 * where \f$i_d=1...N_d\f$ are the indices in each direction \f$d\f$,
 * \f$L_d\f$ the mesh lengths in each direction \f$d\f$ and
 * \f$\{\vec{e}_d,~d=1...d_1\}\f$ an orthonormal basis of
 * \f$\mathbb{R}^{d_1}\f$.
 * 
 * This class also provides a method to flatten a multidimensional index
 * \f$(i_1,i_2,...i_{d_1})\f$ with the following formula:
 * \f[
 * 		\alpha\left[i_1,...i_{d_1}\right] = \sum_{j=1}^{d_1}i_j\,w_j,
 * \f]
 * where the weights \f$w_j\f$ are defined with the following recurrence
 * formula:
 * \f[
 * 		w_j=w_{j+1}\,N_{j+1},\quad w_{d_1}=1
 * \f]
 * See InterpolatedMapping for a concrete use of this class.
 */
template <int dim>
class CartesianMesh {
public:
	/**
	 * Constructor.
	 * 
	 * @param mesh_origin A vector containing the mesh origin
	 * \f$\vec{x}_0\f$.
	 * 
	 * @param mesh_lengths A vector containing the mesh lengths in each
	 * direction \f$d=1...d_1\f$.
	 * 
	 * @param n_points_per_dim A vector containing the numbers of mesh
	 * points in each direction \f$d=1...d_1\f$.
	 */
	CartesianMesh(
		const Vector<dim,double> &mesh_origin,
		const Vector<dim,double> &mesh_lengths,
		const Vector<dim,long> &n_points_per_dim,
		bool last_idx_fastest = true); 

	/**
	 * Returns the flattened index computed with the formula given in
	 * the description of this class.
	 * 
	 * @param indices A vector containing the multi-dimensional index.
	 */
	unsigned long flatten(const Vector<dim,long> &indices) const {
		return (indices,flatten_weights);
	};
	/**
	 * Refine this mesh by dividing by 2 the cell lengths in each
	 * directions.
	 */
	void refine();

	/**
	 * A vector containing the mesh origin \f$\vec{x}_0\f$.
	 */
	Vector<dim,double> origin;
	/**
	 * A vector containing the cell lengths in each direction
	 * \f$d=1...d_1\f$.
	 */
	Vector<dim,double > cell_lengths;
	/**
	 * A vector containing the mesh lengths in each direction
	 * \f$d=1...d_1\f$.
	 */
	Vector<dim,double> lengths;
	/**
	 * A vector containing the weigths \f$w_j\f$ used in the flattening
	 * formula.
	 */
	Vector<dim,long> flatten_weights;
	/**
	 * A vector containing the numbers of mesh points in each direction
	 * \f$d=1...d_1\f$.
	 */
	Vector<dim,long> n_points_per_dim;
	/**
	 * Number of points per cell in dimension <tt>dim_1</tt>.
	 */
	int n_points_per_cell;
	/**
	 * Total number of points in the mesh.
	 */
	long n_points;

	/**
	 * Is the last index of the MultiDimIndex the fastest to change when
	 * we iterate?
	 */
	bool last_idx_fastest;
};

#endif
