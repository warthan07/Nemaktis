#include "ODEFunction.h"

template <int dim>
ODEFunction<dim>::ODEFunction(
		const std::shared_ptr<const DefinitionDomain<dim> > &def_domain) :
	def_domain(def_domain) {}

template <int dim>
ExtraordinaryODEFunction<dim>::ExtraordinaryODEFunction(
		double eps_par, double eps_perp,
		const std::shared_ptr<Mapping<dim,3,double> > &director_field) :
	ODEFunction<dim>(director_field->get_def_domain()),
	eps_par(eps_par),
	eps_perp(eps_perp),
	eps_a(eps_par-eps_perp),
	eps_a_bar( (eps_par*eps_perp) / (eps_par-eps_perp) ),
	director_field(director_field->clone()) {}

template <int dim>
ExtraordinaryODEFunction<dim>::ExtraordinaryODEFunction(
		const ExtraordinaryODEFunction<dim> &ode) :
	ODEFunction<dim>(ode.director_field->get_def_domain()),
	eps_par(ode.eps_par),
	eps_perp(ode.eps_perp),
	eps_a(ode.eps_a),
	eps_a_bar(ode.eps_a_bar),
	director_field(ode.director_field->clone()) {}

template <int dim>
std::shared_ptr<ODEFunction<dim> > ExtraordinaryODEFunction<dim>::clone() const {
	return std::make_shared<ExtraordinaryODEFunction<dim> >(*this);
}

template <int dim>
void ExtraordinaryODEFunction<dim>::compute_ode_function(
		const Ray &ray,
		Vector<3,double> &pos_derivative,
		Vector<3,double> &moment_derivative) {

	// We get the value of the director and its derivative
	Vector<3,double> n;
	if(!director_field->get_value_embedded(ray.cur_pos, n))
		throw(std::string(
			"The ray is outside the definition domain of the "
			"director field (ExtraordinaryODEFunction::compute_ode_function)"));

	Matrix<3,3,double> grad_n;
	director_field->get_gradient_embedded(ray.cur_pos, grad_n);

	// We compute the ode functions
	pos_derivative =
		ray.cur_moment / eps_par
		+ (n,ray.cur_moment) * n / eps_a_bar;
	moment_derivative = 
		- (n,ray.cur_moment) * (grad_n & ray.cur_moment) / eps_a_bar;
}

template <int dim>
void ExtraordinaryODEFunction<dim>::compute_polarisation(
		const Ray &ray,
		Vector<3,std::complex<double> > &pol_res,
		bool forward_search) {

	// We get the current moment 
	Vector<3,std::complex<double> > moment;
	if(ray.is_evanescent)
		moment = ray.complex_moment;
	else
		moment = ray.cur_moment;

	// We get the value of the director
	Vector<3,double> n;
	if(!director_field->get_value_embedded(ray.cur_pos, n))
		throw(std::string(
			"The ray is outside the definition domain of the "
			"director field (ExtraordinaryODEFunction::compute_polarisation)"));

	// We compute the polarisation
	Vector<3,std::complex<double> > tmp = n^moment;
	if(tmp.norm()>1e-10) {
		pol_res =
			eps_perp * n -
			(n,moment) * moment;
		pol_res /= pol_res.norm();
	}
	else if(forward_search) {
		Vector<3,double> x = ray.cur_pos;
		Vector<3,double> old_n = n;
		while(director_field->get_value_embedded(x, n)) {
			tmp = n - (n,old_n) * old_n;
			if(tmp.norm()>1e-10) {
				pol_res = tmp / tmp.norm();
				break;
			}
			else
				x += 0.01*old_n;
		}
	}
}

template <int dim>
void ExtraordinaryODEFunction<dim>::update_optical_index(
		const Ray &ray, double &optical_index_res) {

	// We get the value of the director
	Vector<3,double> n;
	if(!director_field->get_value_embedded(ray.cur_pos, n))
		throw(std::string(
			"The ray is outside the definition domain of the "
			"director field (ExtraordinaryODEFunction::update_optical_index)"));

	// We compute the optical index
	optical_index_res = 1. / std::sqrt(
		1. / eps_par +
		std::pow((n,ray.cur_moment), 2.) / (eps_perp*eps_a_bar));
}

template <int dim>
void ExtraordinaryODEFunction<dim>::update_moment(Ray &ray) {

	Vector<3,double> nu;
	if(!this->def_domain->get_normal_embedded(ray.cur_pos, nu))
		throw(std::string(
			"The ray is outside the definition domain of the "
			"director field"));

	Vector<3,double> n;
	if(!director_field->get_value_embedded(ray.cur_pos, n))
		throw(std::string(
			"The ray is outside the definition domain of the "
			"director field"));

	Vector<3,double> pn = ray.cur_moment - (ray.cur_moment,nu) * nu;
	double n_nu = (n,nu);
	double n_pn = (n,pn);

	double A = eps_perp + eps_a * std::pow(n_nu, 2.);
	double B = eps_a * n_nu * n_pn;
	double C =
		eps_perp * std::pow(pn.norm(), 2.) +
		eps_a * std::pow(n_pn, 2.) - eps_par * eps_perp;

	if(B*B>=A*C) {
		double xi = ( - std::sqrt(B*B-A*C) - B ) / A;
		ray.cur_moment = pn + xi * nu;
	}
	else {
		std::complex<double> xi(-B/A, -std::sqrt(A*C-B*B)/A);
		ray.is_evanescent = true;
		ray.complex_moment = pn + xi * nu;
	}
}

template <int dim>
double ExtraordinaryODEFunction<dim>::hamiltonian_value(Ray &ray) {

	// We get the value of the director
	Vector<3,double> n;
	if(!director_field->get_value_embedded(ray.cur_pos, n)) {
		throw(std::string(
			"The ray is outside the definition domain of the "
			"director field"));
	}

	return ( (ray.cur_moment,ray.cur_moment) / eps_par +
		std::pow((n,ray.cur_moment), 2.) / eps_a_bar ) / 2.;
}

template <int dim>
OrdinaryODEFunction<dim>::OrdinaryODEFunction(
		double eps_perp,
		const std::shared_ptr<Mapping<dim,3,double> > &director_field) :
	ODEFunction<dim>(director_field->get_def_domain()),
	eps_perp(eps_perp),
	director_field(director_field->clone()) {}

template <int dim>
OrdinaryODEFunction<dim>::OrdinaryODEFunction(
		const OrdinaryODEFunction<dim> &ode) :
	ODEFunction<dim>(ode.def_domain),
	eps_perp(ode.eps_perp),
	director_field(ode.director_field->clone()) {}

template <int dim>
std::shared_ptr<ODEFunction<dim> > OrdinaryODEFunction<dim>::clone() const {
	return std::make_shared<OrdinaryODEFunction<dim> >(*this);
}

template <int dim>
void OrdinaryODEFunction<dim>::compute_ode_function(
		const Ray &ray,
		Vector<3,double> &pos_derivative,
		Vector<3,double> &moment_derivative) {

	pos_derivative = ray.cur_moment / eps_perp;
	moment_derivative = 0;
}

template <int dim>
void OrdinaryODEFunction<dim>::compute_polarisation(
		const Ray &ray,
		Vector<3,std::complex<double> > &pol_res,
		bool forward_search) {

	// We get the current moment 
	Vector<3,std::complex<double> > moment;
	if(ray.is_evanescent)
		moment = ray.complex_moment;
	else
		moment = ray.cur_moment;

	// We get the value of the director
	Vector<3,double> n;
	if(!director_field->get_value_embedded(ray.cur_pos, n))
		throw(std::string(
			"The ray is outside the definition domain of the "
			"director field"));

	// We compute the polarisation
	Vector<3,std::complex<double> > tmp = n^moment;
	if(tmp.norm()>1e-10)
		pol_res = tmp / tmp.norm();
	else if(forward_search) {
		Vector<3,double> x = ray.cur_pos;
		Vector<3,double> old_n = n;
		while(director_field->get_value_embedded(x, n)) {
			tmp = n ^ old_n;
			if(tmp.norm()>1e-10) {
				pol_res = tmp / tmp.norm();
				break;
			}
			else
				x += 0.01*old_n;
		}
	}
}

template <int dim>
void OrdinaryODEFunction<dim>::update_moment(Ray &ray) {

	Vector<3,double> nu;
	if(!this->def_domain->get_normal_embedded(ray.cur_pos, nu))
		throw(std::string(
			"The ray is outside the definition domain of the "
			"director field"));

	Vector<3,double> pn = ray.cur_moment - (ray.cur_moment,nu) * nu;
	double pn_2 = std::pow(pn.norm(), 2.);

	if(eps_perp >= pn_2) {
		double xi = -std::sqrt(eps_perp - pn_2);
		ray.cur_moment = pn + xi * nu;
	}
	else {
		std::complex<double> xi(0, -std::sqrt(pn_2 - eps_perp));
		ray.is_evanescent = true;
		ray.complex_moment = pn + xi * nu;
	}
}

template <int dim>
void OrdinaryODEFunction<dim>::update_optical_index(
		const Ray &ray, double &optical_index_res) {

	optical_index_res = std::sqrt(eps_perp);
}

template <int dim>
double OrdinaryODEFunction<dim>::hamiltonian_value(Ray &ray) {

	return (ray.cur_moment,ray.cur_moment) / (2*eps_perp);
}

template <int dim>
IsotropicODEFunction<dim>::IsotropicODEFunction(
		double eps_r,
		const std::shared_ptr<const DefinitionDomain<dim> > &def_domain) :
	ODEFunction<dim>(def_domain),
	eps_r(eps_r) {}

template <int dim>
std::shared_ptr<ODEFunction<dim> > IsotropicODEFunction<dim>::clone() const {
	return std::make_shared<IsotropicODEFunction<dim> >(*this);
}
	
template <int dim>
void IsotropicODEFunction<dim>::compute_ode_function(
		const Ray &ray,
		Vector<3,double> &pos_derivative,
		Vector<3,double> &moment_derivative) {

	pos_derivative = ray.cur_moment / eps_r;
	moment_derivative = 0;
}

template <int dim>
void IsotropicODEFunction<dim>::update_optical_index(
		const Ray &ray, double &optical_index_res) {

	optical_index_res = std::sqrt(eps_r);
}

template <int dim>
void IsotropicODEFunction<dim>::update_moment(Ray &ray) {

	Vector<3,double> nu;
	if(!this->def_domain->get_normal_embedded(ray.cur_pos, nu))
		throw(std::string(
			"You are trying to update the moment of a ray which\n"
			"is not on an interface. Check the definition of the domains."));

	Vector<3,double> pn = ray.cur_moment - (ray.cur_moment,nu) * nu;
	double pn_2 = std::pow(pn.norm(), 2.);

	if(eps_r >= pn_2) {
		double xi = - std::sqrt(eps_r - pn_2);
		ray.cur_moment = pn + xi * nu;
	}
	else {
		std::complex<double> xi(0, -std::sqrt(pn_2 - eps_r));
		ray.is_evanescent = true;
		ray.complex_moment = pn + xi * nu;
	}
}

template <int dim>
double IsotropicODEFunction<dim>::hamiltonian_value(Ray &ray) {

	return (ray.cur_moment,ray.cur_moment) / (2*eps_r);
}

template <int dim>
ODESequence<dim>::ODESequence(std::vector<std::shared_ptr<ODEFunction<dim> > > &_odes) {
	for(auto o : _odes)
		odes.push_back(o->clone());
}

template <int dim>
ODESequence<dim>::ODESequence(const ODESequence<dim> &ode_seq) {
	for(auto o : ode_seq.odes)
		odes.push_back(o->clone());
}

template class ODEFunction<2>;
template class ODEFunction<3>;

template class ExtraordinaryODEFunction<2>;
template class ExtraordinaryODEFunction<3>;

template class OrdinaryODEFunction<2>;
template class OrdinaryODEFunction<3>;

template class IsotropicODEFunction<2>;
template class IsotropicODEFunction<3>;

template class ODESequence<2>;
template class ODESequence<3>;
