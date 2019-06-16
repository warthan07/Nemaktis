#ifndef RAYBUNDLE_H
#define RAYBUNDLE_H

#include <complex>
#include <memory>

#include "linalg.h"

template <int dim>
class RayStepper;

template <int dim>
class FresnelTransmission;

template <int dim>
class ThinLens;

/**
 * A container class for a ray, represented by a point \f$\vec{\eta} = 
 * (\vec{x},\vec{p})\f$ in the phase space, where \f$\vec{x}\f$ is the
 * current position and \f$\vec{p}\f$ is the conjugate moment.
 */
class Ray {
public:
	/**
	 * Constructor.
	 * 
	 * @param init_pos Vector containing the starting point of the ray.
	 * 
	 * @param init_moment Vector containing the starting moment of the
	 * ray.
	 * 
	 * @param init_optical_index Value of the optical index at the
	 * starting point.
	 */
	Ray(
		const Vector<3,double> &init_pos, const Vector<3,double> &init_moment);

	/**
	 * The current spatial position.
	 */
	Vector<3,double> cur_pos;
	/**
	 * The current moment.
	 */
	Vector<3,double> cur_moment;
	/**
	 * The optical length
	 * \f$\bar{s}\equiv\int_0^s\vec{p}\cdot\mathrm{d}\vec{r}\f$ between
	 * the initial shooting point and current point. For an isotropic
	 * medium of optical index \f$n\f$, \f$\bar{s}\equiv n\,s\f$.
	 */
	double cur_optical_length;

	/**
	 * A boolean indicating if this ray has become evanescent after a
	 * reflection or transmission at an interface.
	 */
	bool is_evanescent;
	/**
	 * If this ray is evanescent, the complex wave moment is stored in
	 * this attribute
	 */
	Vector<3,std::complex<double> > complex_moment;
};

/**
 * A bundle of (dim+1) rays allowing to easily compute a finite-difference
 * approximation of the geometrical spreading \f$q\f$, defined with the
 * following formula:
 * \f[
 * 		q = \det\left[\frac{\partial \vec{x}}{\partial\vec{x}^{(i)}}
 * 			\left(\bar{s},\vec{x}^{(i)}\right)\right]
 * \f]
 * where \f$\vec{x}^{(i)}\f$ is the starting point of the ray, \f$\bar{s}\f$ is the
 * current optical length of the ray and \f$\vec{x}\f$ is the current
 * end point of the ray.
 * 
 * To this end, this class allows to store the trajectories
 * of a central rays \f$\vec{x}_0=\vec{x}_f(\bar{s},\vec{x}^{(i)})\f$ and
 * \f$\rm dim\f$ satellite rays
 * \f$\vec{x}_k=\vec{x}_f(\bar{s},\vec{x}^{(i)}+\mu\,\vec{e}_k)\f$,
 * where \f$\left\{\vec{e}_k,\,k=1\,...\,\rm dim\right\}\f$ is an
 * orthonormal basis for the physical space. The finite difference step
 * \f$\mu\f$ is stored in \link fd_step \endlink, while the rays
 * trajectories are stored in the array \link rays \endlink.
 */
template <int dim>
class RayBundle {
public:
	/**
	 * Constructor.
	 * 
	 * @param init_pos Vector containing the starting point of the ray.
	 * 
	 * @param init_moment Vector containing the starting moment of the
	 * ray.
	 * 
	 * @param fd_step Finite-difference step \f$\mu\f$ used to generate the
	 * satellite rays from the central ray.
	 * 
	 * @param wavelength The wavelength \f$\lambda\f$ of the ray.
	 * 
	 * @param init_amplitude Amplitude \f$E\f$ at the starting point.
	 * 
	 * @param init_optical_index Optical index \f$n_{\rm eff}\f$ at the
	 * starting point (see \link Ray::cur_optical_index \endlink for a
	 * precise definition of \f$n_{\rm eff}\f$).
	 */
	RayBundle(
		const Vector<3,double> &init_pos, const Vector<3,double> &init_moment,
		const Vector<3,double> (&init_pol)[2],	double init_optical_index,
		double init_ampl, double fd_step);

	void project_on_analyzer(Vector<3,double> &analyzer);

	/**
	 * Returns the current value of the electric field amplitude, computed
	 * with the following formula:
	 * \f[
	 * 		E=\frac{T\,J^{(i)}_{\rm eff}}{n_{\rm eff}\sqrt{|q|}},
	 * \f]
	 * where \f$J^{(i)}_E=n^{(i)}_{\rm eff}\sqrt{q^{(i)}}E^{(i)}\f$ is the 
	 * rescaled amplitude at the starting point (\link
	 * #init_rescaled_amplitude \endlink), \f$T\f$ is the product of
	 * all transmission factor at the crossed interfaces (\link
	 * #cur_transmission_factor \endlink), \f$n_{\rm eff}\f$ is the
	 * current optical index (\link Ray::cur_optical_index \endlink),
	 * and \f$q\f$ is the current geometrical spreading (\link
	 * #cur_geom_spreading \endlink).
	 * 
	 * This formula comes the fact that \f$J_E=n_{\rm eff}\sqrt{q}E\f$ is
	 * conserved along a ray if there is no discontinuity of the optical
	 * index. In presence of discontinuities, \f$J_E\f$ is only constant
	 * by part, and transmission factors must be computed with the
	 * Fresnel formula.
	 *
	 * Returns the current value of the electric field phase, computed
	 * with the following formula:
	 * \f[
	 * 		\phi = \frac{2\pi\,\bar{s}}{\lambda}+\frac{m\,\pi}{2}.
	 * \f]
	 * Here, \f$\bar{s}\f$ is the current optical length of the ray
	 * (\link Ray::cur_optical_length \endlink), \f$\lambda\f$ is the
	 * wavelength (\link wavelength \endlink), and \f$m\f$ is the number of
	 * times the ray has crossed a caustic (\link cur_caustic_number
	 * \endlink).
	 */
	void compute_E_field(Vector<3,std::complex<double> > &res, int pol_idx);
	void compute_B_field(Vector<3,std::complex<double> > &res, int pol_idx);
	double optical_length();

	bool is_NA_limited();
	bool is_evanescent();

	/**
	 * Returns the current end point of the i-th ray.
	 */
	Vector<3,double>& cur_pos(int i=0) {
		return rays[i]->cur_pos;
	}
	/**
	 * Returns the current moment of the i-th ray.
	 */
	Vector<3,double>& cur_moment(int i=0) {
		return rays[i]->cur_moment;
	}
	double get_fd_step() {
		return fd_step;
	}

	friend class RayStepper<dim>;
	friend class FresnelTransmission<dim>;
	friend class ThinLens<dim>;

	/**
	 * Update the value of \link cur_geom_spreading \endlink with the
	 * current value of the geometrical spreading, using the
	 * following finite difference approximation:
	 * \f[
	 * 		q = \det{\left(
	 * 			\frac{\vec{x}_1-\vec{x}_0}{\mu}\,\Bigg|\,
	 * 			\frac{\vec{x}_2-\vec{x}_0}{\mu}\,\Bigg|\,
	 * 			\frac{\vec{x}_3-\vec{x}_0}{\mu}
	 * 		\right)},
	 * \f]
	 * where \f$\vec{x}_0=\vec{x}_f(\bar{s},\vec{x}^{(i)})\f$,
	 * \f$\vec{x}_k=\vec{x}_f(\bar{s},\vec{x}^{(i)}+\mu\,\vec{e}_k)\f$,
	 * and \f$\left\{\vec{e}_k,\,k=1\,...\,\rm dim\right\}\f$
	 * is an orthonormal basis for the physical space.
	 * 
	 * In addition, this method test if \f$q\f$ changed its sign (i.e.
	 * the ray has crossed a caustic) since the last update, in which
	 * case \link #cur_caustic_number \endlink is incremented.  This
	 * method is only approximate since it is possible that the
	 * geometrical spreading changes its sign two times during one euler
	 * update. For an exact update of cur_caustic_number, use the method
	 * \link exact_geometrical_spreading_update \endlink.
	 */
	void simple_geometrical_spreading_update(
		bool backtrack = false);
	void exact_geometrical_spreading_update(
		Vector<3,double> (&pos_derivatives)[dim+1], double step_length);

private:
	double compute_geometrical_spreading();

	/**
	 * Array containing pointers to the central ray and 4 satellite rays.
	 */
	std::shared_ptr<Ray> rays[dim+1];

	/**
	 * The current geometrical spreading \f$q\f$.
	 */
	double cur_geom_spreading;
	/**
	 * Current value of the jacobian matrix of the transformation
	 * \f$\vec{x}^{(i)}\rightarrow\vec{x}(\bar{s},\vec{x}^{(i)})\f$.
	 */
	Matrix<dim,dim,double> cur_jac;
	/**
	 * The current transmission factors \f$T\f$ for each vectors of the
	 * polarisation basis.
	 */
	std::complex<double> cur_transmission_factors[2];
	/**
	 * The current optical index \f$n_{\rm eff}\f$, which depends on the
	 * type of ray.  
	 * 
	 * In an isotropic medium, \f$n_{\rm
	 * eff}=\sqrt{\epsilon_r}\f$. For an ordinary ray, \f$n_{\rm
	 * eff}=\sqrt{\epsilon_\perp}\f$. For an extraordinary ray, 
	 * \f$n_{\rm eff}=\sqrt{\epsilon_\parallel-
	 * \epsilon_a\left(\vec{n}\cdot\vec{u}_s\right)^2}\f$, where 
	 * \f$\vec{n}\f$ is the director and \f$\vec{u}_s\f$ is the normed ray
	 * tangent vector. 
	 * 
	 * In all cases, the effective optical index can be
	 * computed with the following formula: \f$n_{\rm
	 * eff}=\frac{\mathrm{d}\bar{s}}{\mathrm{d}s}=\vec{p}\cdot\vec{u}_s\f$.
	 */
	double cur_optical_index;
	/** 
	 * The current polarisation basis.
	 */
	Vector<3,std::complex<double> > cur_pols[2];

	/**
	 * Is the current ray bundle limited by the numerical aperture of
	 * the lens?
	 */
	bool NA_limited;
	/**
	 * Is the current ray bundle became evanescent at some point?
	 */
	bool evanescent;

	/**
	 * The finite-difference step \f$\mu\f$ used to compute the geometrical
	 * spreading. 
	 */
	const double fd_step;
	/**
	 * The rescaled amplitude \f$J_E=n_{\rm eff}\sqrt{q}E\f$
	 * at the starting point.
	 */
	const double init_rescaled_amplitude;

};

#endif
