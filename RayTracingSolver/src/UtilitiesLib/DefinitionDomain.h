#ifndef DOMAININTERFACE_h
#define DOMAININTERFACE_h

#define SURFACE_EPS 1e-12

#include <memory>

#include "linalg.h"
#include "error.h"
#include "CartesianMesh.h"

/**
 * Pure virtual class which specifies a definition domain
 * \f$\Omega\subset\mathbb{R}^d\f$ for a Mapping object.
 */
template <int dim>
class DefinitionDomain {
public:
	/**
	 * Constructor.
	 */
	DefinitionDomain() {}; 

	/**
	 * Return true (false) if <tt>pos</tt> is inside (outside) the
	 * domain specified by this class.
	 *
	 * @param pos The spatial position to be tested.
	 */
	virtual bool is_inside(
		const Vector<dim,double> &pos, double eps = SURFACE_EPS) const = 0;
	virtual bool is_inside_embedded(
		const Vector<3,double> &pos, double eps = SURFACE_EPS) const;

	/**
	 * Return true (false) if <tt>pos</tt> is on the interface
	 * delimiting \f$\Omega\f$ within an absolute precision of
	 * <tt>eps</tt>
	 *
	 * @param pos The spatial position to be tested.
	 * 
	 * @param eps The absolute precision for this test.
	 */
	virtual bool is_on_surface(
		const Vector<dim,double> &pos, double eps = SURFACE_EPS) const = 0;
	virtual bool is_on_surface_embedded(
		const Vector<3,double> &pos, double eps = SURFACE_EPS) const;
	/**
	 * Compute the interface normal in <tt>pos</tt> and give the
	 * result by reference in <tt>res</tt>. If <tt>pos</tt> is not on
	 * the interface, false is returned, else true is returned.
	 *
	 * @param pos The spatial position at which the normal must be
	 * computed.
	 *
	 * @param res Reference to the result vector which will store the
	 * value of the normal.
	 * 
	 * @param eps The absolute precision for the test \link
	 * is_on_surface \endlink.
	 */
	virtual bool get_normal(
		const Vector<dim,double> &pos, Vector<dim,double> &res,
		double eps = SURFACE_EPS) const = 0;
	virtual bool get_normal_embedded(
		const Vector<3,double> &pos, Vector<3,double> &res,
		double eps = SURFACE_EPS) const;
	/**
	 * Compute the nearest intersection between a straight line and the
	 * interface and give the result by reference in <tt>res</tt>. The
	 * straight line is specified by an initial position <tt>pos</tt>
	 * and a velocity vector <tt>vel</tt>. If an intersection is found,
	 * this function returns the time of travel between the starting
	 * point and the intersection point, else it returns -1.
	 *
	 * @param pos Starting point of the straight line.
	 *
	 * @param vel Velocity on the straight line.
	 * 
	 * @param res Reference to the result vector which will store the
	 * value of the nearest intersection (if found).
	 */
	virtual double compute_intersection(
		const Vector<dim,double> &pos, const Vector<dim,double> &vel,
		Vector<dim,double> &res) const = 0;
	virtual double compute_intersection_embedded(
		const Vector<3,double> &pos, const Vector<3,double> &vel,
		Vector<3,double> &res) const;
	/**
	 * Compute the nearest orthogonal projection of a point on the domai
	 * interface, give the result by reference in <tt>res</tt>, and
	 * return the distance between the given point and its projection.
	 *
	 * @param pos Input point.
	 *
	 * @param res Reference to the result vector which will store the
	 * value of the nearest projection on the interface.
	 */
	virtual double compute_projection(
		const Vector<dim,double> &pos, Vector<dim,double> &res,
		double eps = 0) const = 0;
	virtual double compute_projection_embedded(
		const Vector<3,double> &pos, Vector<3,double> &res,
		double eps = 0) const;
};

/**
 * A specialization of DefinitionDomain when the domain \f$\Omega\f$ is a
 * ball, specified by its center \f$\vec{x}_0\f$ and its radius
 * \f$R\f$.
 */
template <int dim>
class SphericalDomain : public DefinitionDomain<dim> {
public:
	/**
	 * Constructor.
	 * 
	 * @param sphere_center The sphere center \f$\vec{x}_0\f$.
	 * 
	 * @param sphere_radius The sphere radius \f$R\f$.
	 */
	SphericalDomain(
		const Vector<dim,double> &sphere_center,
		double sphere_radius);

	virtual bool is_inside(
		const Vector<dim,double> &pos, double eps = SURFACE_EPS) const;
	virtual bool is_on_surface(
		const Vector<dim,double> &pos, double eps = SURFACE_EPS) const;
	virtual bool get_normal(
		const Vector<dim,double> &pos, Vector<dim,double> &res,
		double eps = SURFACE_EPS) const;
	virtual double compute_intersection(
		const Vector<dim,double> &pos, const Vector<dim,double> &vel,
		Vector<dim,double> &res) const;
	virtual double compute_projection(
		const Vector<dim,double> &pos, Vector<dim,double> &res,
		double eps = 0) const;

private:
	/**
	 * The sphere center \f$\vec{x}_0\f$.
	 */
	const Vector<dim,double> center;
	/**
	 * The sphere radius \f$R\f$.
	 */
	const double radius;
};

/**
 * A specialization of DefinitionDomain when the domain \f$\Omega\f$ is a
 * half-disk, specified by its radius \f$R\f$. The straight side of the
 * half-disk is aligned with the y-axis, and is vertically aligned.
 */
class HalfDiskDomain : public DefinitionDomain<2> {
public:
	/**
	 * Constructor.
	 * 
	 * @param disk_radius The half-disk radius \f$R\f$.
	 */
	HalfDiskDomain(double disk_radius);

	virtual bool is_inside(
		const Vector<2,double> &pos, double eps = SURFACE_EPS) const;
	virtual bool is_on_surface(
		const Vector<2,double> &pos, double eps = SURFACE_EPS) const;
	virtual bool get_normal(
		const Vector<2,double> &pos, Vector<2,double> &res,
		double eps = SURFACE_EPS) const;
	virtual double compute_intersection(
		const Vector<2,double> &pos, const Vector<2,double> &vel,
		Vector<2,double> &res) const;
	virtual double compute_projection(
		const Vector<2,double> &pos, Vector<2,double> &res,
		double eps = 0) const;

private:
	/**
	 * The disk radius \f$R\f$.
	 */
	const double radius;
};

/**
 * A specialization of DefinitionDomain when the domain \f$\Omega\f$ is a
 * rectangular parallelotope specified with the following points:
 * \f[
 * 		\vec{x}\left[i_1,...i_{d_1}\right] = \vec{x}_0 + \sum_{d=1}^{d_1}
 * 			\frac{i_d-1}{N_d-1}\,L_d\,\vec{e}_d,
 * \f]
 * where \f$i_d=1...2\f$ are the point indices in each direction \f$d\f$,
 * \f$L_d\f$ the mesh lengths in each direction \f$d\f$ and
 * \f$\{\vec{e}_d,~d=1...d_1\}\f$ an orthonormal basis of
 * \f$\mathbb{R}^{d_1}\f$, with \f$d_1=\f$<tt>dim</tt>.
 */
template <int dim>
class ParallelotopeDomain : public DefinitionDomain<dim> {
public:
	/**
	 * Default constructor.
	 * 
	 * @param domain_origin The parallelotope origin \f$\vec{x}_0\f$.
	 * 
	 * @param domain_lengths The parallelotope lengths \f$L_d\f$ in each
	 * direction \f$d=1...d_1\f$.
	 */
	ParallelotopeDomain(
		const Vector<dim,double> &domain_origin,
		const Vector<dim,double> &domain_lengths);
	/**
	 * Constructor allowing to get the same definition domain as a
	 * CartesianMesh object.
	 * 
	 * @param mesh The CartesianMesh object which will initialize this
	 * class.
	 */
	ParallelotopeDomain(
		const CartesianMesh<dim> &mesh);


	virtual bool is_inside(
		const Vector<dim,double> &pos, double eps = SURFACE_EPS) const;
	virtual bool is_on_surface(
		const Vector<dim,double> &pos, double eps = SURFACE_EPS) const;
	virtual bool get_normal(
		const Vector<dim,double> &pos, Vector<dim,double> &res,
		double eps = SURFACE_EPS) const;
	virtual double compute_intersection(
		const Vector<dim,double> &pos, const Vector<dim,double> &vel,
		Vector<dim,double> &res) const;
	virtual double compute_projection(
		const Vector<dim,double> &pos, Vector<dim,double> &res,
		double eps = 0) const;

private:
	/**
	 * A vector containing the domain origin \f$\vec{x}_0\f$.
	 */
	Vector<dim,double> domain_origin;
	/**
	 * A vector containing the domain lengths \f$L_d\f$ in each
	 * direction \f$d=1...d_1\f$.
	 */
	Vector<dim,double> domain_lengths;

};

/**
 * A specialization of DefinitionDomain when the domain \f$\Omega\f$ is a
 * slab delimited by the two following hyper-planes:
 * \f[
 * 		\vec{x}\cdot\vec{\nu} = h_1,
 * \f]
 * \f[
 * 		\vec{x}\cdot\vec{\nu} = h_2,
 * \f]
 * where \f$\vec{\nu}\f$ is a unit vector orthogonal to the two
 * hyper-planes and \f$h_1<h_2\f$.
 * 
 */
template <int dim>
class SlabDomain : public DefinitionDomain<dim> {
public:
	/**
	 * Constructor.
	 * 
	 * Note that this constructor will swap the definition of
	 * <tt>h_1</tt> and <tt>h_2</tt> if <tt>h_1>h_2</tt>, in order to
	 * enforce the condition \f$h_1<h_2\f$.
	 * 
	 * @param slab_normal The slab normal \f$\vec{\nu}\f$.
	 * 
	 * @param h_1 The parameter \f$h_1\f$ of the first hypersurface
	 * delimiting the slab.
	 * 
	 * @param h_2 The parameter \f$h_2\f$ of the second hypersurface
	 * delimiting the slab.
	 */
	SlabDomain(const Vector<dim,double> &slab_normal, double h_1, double h_2);

	virtual bool is_inside(
		const Vector<dim,double> &pos, double eps = SURFACE_EPS) const;
	virtual bool is_on_surface(
		const Vector<dim,double> &pos, double eps = SURFACE_EPS) const;
	virtual bool get_normal(
		const Vector<dim,double> &pos, Vector<dim,double> &res,
		double eps = SURFACE_EPS) const;
	virtual double compute_intersection(
		const Vector<dim,double> &pos, const Vector<dim,double> &vel,
		Vector<dim,double> &res) const;
	virtual double compute_projection(
		const Vector<dim,double> &pos, Vector<dim,double> &res,
		double eps = 0) const;

private:
	/**
	 * A vector containing the slab normal \f$\vec{\nu}\f$.
	 */
	Vector<dim,double> slab_normal;
	/**
	 * The parameter \f$h_1\f$ of the first hypersurface delimiting the
	 * slab.
	 */
	double h_1;
	/**
	 * The parameter \f$h_2\f$ of the second hypersurface delimiting the
	 * slab.
	 */
	double h_2;
};

/**
 * A specialization of DefinitionDomain which allows to substract two
 * arbitrary definition domain.
 */
template <int dim>
class SubstractedDomain : public DefinitionDomain<dim> {
public:
	/**
	 * Constructor.
	 */
	SubstractedDomain(
		const std::shared_ptr<DefinitionDomain<dim> > &outer_domain,
		const std::shared_ptr<DefinitionDomain<dim> > &inner_domain);

	virtual bool is_inside(
		const Vector<dim,double> &pos, double eps = SURFACE_EPS) const;
	virtual bool is_on_surface(
		const Vector<dim,double> &pos, double eps = SURFACE_EPS) const;
	virtual bool get_normal(
		const Vector<dim,double> &pos, Vector<dim,double> &res,
		double eps = SURFACE_EPS) const;
	virtual double compute_intersection(
		const Vector<dim,double> &pos, const Vector<dim,double> &vel,
		Vector<dim,double> &res) const;
	virtual double compute_projection(
		const Vector<dim,double> &pos, Vector<dim,double> &res,
		double eps = 0) const;

private:
	/**
	 * Pointer to the outer domain.
	 */
	std::shared_ptr<DefinitionDomain<dim> > outer_domain;
	/**
	 * Pointer to the inner domain.
	 */
	std::shared_ptr<DefinitionDomain<dim> > inner_domain;

};

#endif
