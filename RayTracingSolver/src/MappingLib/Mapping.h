#ifndef MAPPING_H
#define MAPPING_H

#include <memory>
#include <vector>
#include <string>

#include "linalg.h"
#include "CartesianMesh.h"
#include "DefinitionDomain.h"

/**
 * Pure virtual class for a generic mapping
 * \f$\mathbb{R}^{d_1}\rightarrow T^{d_2}\f$, where \f$T=\mathbb{C}\f$
 * or \f$\mathbb{R}\f$ and \f$d_i=\f$<tt>dim_i</tt>. The domain of
 * definition of this mapping is defined with a DomainDefinition object.
 */
template <int dim1, int dim2, typename T>
class Mapping {
public:
	/**
	 * Constructor with a custom DefinitionDomain object.
	 * 
	 * @param def_domain Smart pointer to the DefinitionDomain
	 * object specifying the domain of definition of this mapping.
	 */
	Mapping(const std::shared_ptr<DefinitionDomain<dim1> > &def_domain);

	/**
	 * Clone method, which returns a shared_ptr pointing to a copy of
	 * this object.
	 */
	virtual std::shared_ptr<Mapping<dim1,dim2,T> > clone() const = 0;
	
	/**
	 * Pure virtual method, whose goal is to compute the value
	 * of the interpolated function at a given point <tt>x</tt>.
	 * 
	 * @param x Constant reference to the interpolation point.
	 * 
	 * @param res Reference to the result vector which will store the
	 * interpolated value of the mapping.
	 */
	virtual bool get_value(
		const Vector<dim1,double> &x, Vector<dim2,T> &res) = 0;
	virtual bool get_value_embedded(
		const Vector<3,double> &x, Vector<dim2,T> &res);
	/**
	 * Pure virtual method, whose goal is to compute the gradient
	 * of the interpolated function at a given point <tt>x</tt>.
	 * 
	 * @param x Constant reference to the interpolation point.
	 * 
	 * @param res Reference to the result matrix which will store the
	 * interpolated gradient of the mapping.
	 */
	virtual bool get_gradient(
		const Vector<dim1,double> &x, Matrix<dim1,dim2,T> &res) = 0;
	virtual bool get_gradient_embedded(
		const Vector<3,double> &x, Matrix<3,dim2,T> &res);
	/**
	 * Returns a const reference of the DefinitionDomain smart pointer
	 * def_domain.
	 */
	const std::shared_ptr<DefinitionDomain<dim1> >& get_def_domain() {
		return def_domain;
	}

protected:
	/**
	 * Smart pointer to the DefinitionDomain object specifying the
	 * definition domain of this mapping.
	 */
	std::shared_ptr<DefinitionDomain<dim1> > def_domain;
};

/**
 * Implementation of a mapping \f$\mathbb{R}^{d_1}\rightarrow T^{d_2}\f$
 * (with \f$d_i=\f$<tt>dim_i</tt> and \f$T=\mathbb{R}\f$ or
 * \f$\mathbb{C}\f$) based on lagrange polynomial interpolation for the
 * function values and gradients. This class is purely virtual, the
 * assembly of the interpolation kernel, which depends on the
 * interpolation order, is done on derived classes.
 * 
 * The multidimensional data used for the interpolation is assumed to
 * correspond to the function values
 * \f$\vec{f}\left(\vec{x}\left[i_1,...i_{d_1}\right]\right)\in
 * \mathbb{T}^{d_2}\f$ on the nodes
 * \f$\vec{x}\left[i_1,...i_{d_1}\right]\f$ of a cartesian mesh, with :
 * \f[
 * 		\vec{x}\left[i_1,...i_{d_1}\right] = \vec{x}_0 + \sum_{d=1}^{d_1}
 * 			\frac{i_d-1}{N_d-1}\,L_d\,\vec{e}_d,
 * \f]
 * \f$i_d=1...N_d\f$ the indices in each direction \f$d\f$,
 * \f$L_d\f$ the mesh lengths in each direction \f$d\f$ and
 * \f$\{\vec{e}_d,~d=1...d_1\}\f$ an orthonormal basis of
 * \f$\mathbb{R}^{d_1}\f$.
 * 
 * Data must be provided as a 1D vector of function values:
 * \f[
 * 		\left\{f\left(\vec{x}[\alpha]\right),\quad
 * 		\alpha=1\;...\left(\prod_{d=1}^{d_1}N_d\right)\right\},
 * \f]
 * where we flattened the multidimensional index
 * \f$(i_1,i_2,...i_{d_1})\f$ with the following formula:
 * \f[
 * 		\alpha\left[i_1,...i_{d_1}\right] = \sum_{j=1}^{d_1}i_j\,w_j
 * \f]
 * The weights \f$w_j\f$ are defined with the following recurrence
 * formula:
 * \f[
 * 		w_j=w_{j+1}\,N_{j+1},\quad w_{d_1}=1
 * \f]
 * This choice of flattened index provides the fastest access when
 * iterating on the last index \f$i_{d_1}\f$, which by convention
 * corresponds to the main propagation axis for the light. 
 *
 * Finally, the user has the possibility to restrict the domain of
 * definition of a mapping by giving to this class a DefinitionDomain
 * object (with any constructor of this class). In this case, data is
 * extrapolated near the boundary of the DefinitionDomain in order to
 * provide accurate interpolation even on cells overlapping the dopmain
 * interface.
 */
template <int dim1, int dim2, typename T>
class InterpolatedMapping : public Mapping<dim1,dim2,T> {
public:
	/**
	 * Default constructor.
	 * 
	 * Allows to directly input the data with a vector.
	 * 
	 * @param values Smart pointer to the vector of function values
	 * \f$\vec{f}_\alpha\f$.
	 *
	 * @param mesh SmartPointer to the CartesianMesh object specifying
	 * the mesh on which the interpolation is done.
	 * 
	 * @param def_domain Smart pointer to the DefinitionDomain
	 * object specifying the domain of definition of this mapping.
	 *
	 * @param pol_order The order of interpolation.
	 */
	InterpolatedMapping(
		const std::shared_ptr<std::vector<Vector<dim2,T> > > &values,
		const std::shared_ptr<CartesianMesh<dim1> > &mesh,
		const std::shared_ptr<DefinitionDomain<dim1> > &def_domain,
		int pol_order);

	/**
	 * Compute the value of the interpolated function at a given point
	 * <tt>x</tt> using the values of the interpolating polynome
	 * computed with <tt>assemble_pol_weights</tt>.
	 * 
	 * @param x Constant reference to the interpolation point.
	 * 
	 * @param res Reference to the result vector which will store the
	 * interpolated value of the mapping.
	 */
	virtual bool get_value(
		const Vector<dim1,double> &x, Vector<dim2,T> &res);
	/**
	 * Compute the gradient of the interpolated function at a given
	 * point <tt>x</tt> using the values of the interpolating polynome
	 * computed with <tt>assemble_pol_weights</tt>.
	 *
	 * @param x Constant reference to the interpolation point.
	 * 
	 * @param res Reference to the result matrix which will store the
	 * interpolated gradient of the mapping.
	 */
	virtual bool get_gradient(
		const Vector<dim1,double> &x, Matrix<dim1,dim2,T> &res);
	/**
	 * Normalize the values of the mapping.
	 */
	void normalize();

	/**
	 * Extrapolate the data on each points associated with the indices 
	 * in \link nearby_point_indices \endlik.
	 */
	void extrapolate_data();


protected:
	/**
	 * This method fill the array \link nearby_point_indices \endlink
	 * with the flattened indices of all external points which need
	 * extrapolated data in order to get transparent interpolation near
	 * the domain interface. These points are found by selecting all
	 * points external to the domain but contained in a N by N array of
	 * cells centered on any boundary cells overlapping with the domain
	 * interface, where N is the order of the interpolating polynomial
	 * (1 for linear interpolation, 3 for cubic interpolation).
	 */
	void get_nearby_points_indices();

	/**
	 * A private virtual method which assemble the interpolating
	 * polynomial weights for a given cell. Implemented in child classes
	 * depending on the interpolation order.
	 */
	virtual void assemble_pol_weights(MultiDimIndex<dim1> &cell_origin_indices) = 0;

	/**
	 * An array containing the flattened indices of all external points
	 * near the boundrary which needs extrapolated data, together with
	 * the indices of the nearest internal points for each external
	 * point.
	 */
	std::vector<std::pair<MultiDimIndex<dim1>,MultiDimIndex<dim1> > >
		nearby_point_data;

	/**
	 * The interpolating polynomial order, which will be set when the
	 * constructor of this class is called from a child implementation.
	 */
	int pol_order;

	/**
	 * Smart pointer to the flattened array containing the function
	 * values on each mesh nodes.
	 */
	std::shared_ptr<std::vector<Vector<dim2,T> > > values;


	/**
	 * Flattened index of the last interpolation cell origin.
	 */
	long last_cell_origin_idx;
	/**
	 * An array of vectors used internally to compute the polynomial
	 * weights.
	 */
	std::vector<Vector<dim2,T> > v;
	/**
	 * An array containing the interpoalting polynomial weights for the
	 * current cell. Its size depends on the interpolation order and the
	 * dimension.
	 */
	std::vector<Vector<dim2,T> > pol_coefs;

	/**
	 * A smart pointer to the CartesianMesh object on which the
	 * interpolation is done.
	 */
	std::shared_ptr<const CartesianMesh<dim1> > mesh;
};

#endif
