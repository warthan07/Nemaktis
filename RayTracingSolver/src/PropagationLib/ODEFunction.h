#ifndef ODEFUNCTION_H
#define ODEFUNCTION_H

#include <memory>

#include "RayBundle.h"
#include "Mapping.h"

template <int dim>
class ODEFunction {
public:
	ODEFunction(const std::shared_ptr<const DefinitionDomain<dim> > &def_domain);

	virtual std::shared_ptr<ODEFunction<dim> > clone() const = 0;

	virtual void compute_ode_function(
			const Ray &ray,
			Vector<3,double> &pos_derivative,
			Vector<3,double> &moment_derivative) = 0;
	virtual void compute_polarisation(
			const Ray &ray,
			Vector<3,std::complex<double> > &pol_res,
			bool forward_search = false) = 0;
	virtual void update_optical_index(
			const Ray &ray,
			double &optical_index_res) = 0;
	virtual void update_moment(Ray &ray) = 0;
	virtual double hamiltonian_value(Ray &ray) = 0;
	virtual bool is_isotropic() const = 0;

	std::shared_ptr<const DefinitionDomain<dim> > def_domain;
};

template <int dim>
class ExtraordinaryODEFunction : public ODEFunction<dim> {
public:
	ExtraordinaryODEFunction(
		double eps_par, double eps_perp,
		const std::shared_ptr<Mapping<dim,3,double> > &director_field);
	ExtraordinaryODEFunction(
		const ExtraordinaryODEFunction<dim> &ode);

	virtual std::shared_ptr<ODEFunction<dim> > clone() const;

	virtual void compute_ode_function(
			const Ray &ray,
			Vector<3,double> &pos_derivative,
			Vector<3,double> &moment_derivative);
	virtual void compute_polarisation(
			const Ray &ray,
			Vector<3,std::complex<double> > &pol_res,
			bool forward_search = false);
	virtual void update_optical_index(
			const Ray &ray,
			double &optical_index_res);
	virtual void update_moment(Ray &ray);
	virtual double hamiltonian_value(Ray &ray);
	virtual bool is_isotropic() const { return false; }


private:
	double eps_par;
	double eps_perp;
	double eps_a;
	double eps_a_bar;

	std::shared_ptr<Mapping<dim,3,double> > director_field;
};

template <int dim>
class OrdinaryODEFunction : public ODEFunction<dim> {
public:
	OrdinaryODEFunction(
		double eps_perp,
		const std::shared_ptr<Mapping<dim,3,double> > &director_field);
	OrdinaryODEFunction(
		const OrdinaryODEFunction<dim> &ode);

	virtual std::shared_ptr<ODEFunction<dim> > clone() const;

	virtual void compute_ode_function(
			const Ray &ray,
			Vector<3,double> &pos_derivative,
			Vector<3,double> &moment_derivative);
	virtual void compute_polarisation(
			const Ray &ray,
			Vector<3,std::complex<double> > &pol_res,
			bool forward_search = false);
	virtual void update_optical_index(
			const Ray &ray,
			double &optical_index_res);
	virtual void update_moment(Ray &ray);
	virtual double hamiltonian_value(Ray &ray);
	virtual bool is_isotropic() const { return false; }

private:
	double eps_perp;
	std::shared_ptr<Mapping<dim,3,double> > director_field;
};

template <int dim>
class IsotropicODEFunction : public ODEFunction<dim> {
public:
	IsotropicODEFunction(
		double eps_r,
		const std::shared_ptr<const DefinitionDomain<dim> > &def_domain);

	virtual std::shared_ptr<ODEFunction<dim> > clone() const;

	virtual void compute_ode_function(
			const Ray &ray,
			Vector<3,double> &pos_derivative,
			Vector<3,double> &moment_derivative);
	virtual void compute_polarisation(
			const Ray &ray,
			Vector<3,std::complex<double> > &pol_res,
			bool forward_search = false) {}
	virtual void update_optical_index(
			const Ray &ray,
			double &optical_index_res);
	virtual void update_moment(Ray &ray);
	virtual double hamiltonian_value(Ray &ray);
	virtual bool is_isotropic() const { return true; }

private:
	double eps_r;
};

/**
 * A class representing a sequence of ODEs associated with each medium
 * in which a ray family will travel. Allow deep copy for OpenMP.
 */
template <int dim>
class ODESequence {
public:
	ODESequence() {}
	ODESequence(std::vector<std::shared_ptr<ODEFunction<dim> > > &_odes);
	ODESequence(const ODESequence<dim> &ode_seq);

	std::vector<std::shared_ptr<ODEFunction<dim> > > odes;
};

#endif
