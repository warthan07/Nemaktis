#include <cmath>
#include <iostream>
#include <string>

#include "RayBundle.h"
#include "error.h"

Ray::Ray(
		const Vector<3,double> &init_pos, const Vector<3,double> &init_moment) :
	cur_pos(init_pos),
	cur_moment(init_moment),
	cur_optical_length(0),
	is_evanescent(false) {}

template <int dim>
RayBundle<dim>::RayBundle(
		const Vector<3,double> &init_pos, const Vector<3,double> &init_moment,
		const Vector<3,double> (&init_pols)[2], double init_optical_index,
		double init_ampl, double fd_step) :
	cur_geom_spreading(1),
	cur_optical_index(init_optical_index),
	fd_step(fd_step),
	init_rescaled_amplitude(init_ampl*init_optical_index),
	NA_limited(false),
	evanescent(false) {

	for(int i=0; i<2; i++) {
		cur_pols[i] = init_pols[i];
		cur_transmission_factors[i] = 1;
	}

	rays[0] = std::make_shared<Ray>(init_pos, init_moment);

	Vector<3,double> init_pos_i;
	for(unsigned int i=0; i<dim; ++i) {
		init_pos_i = init_pos;
		init_pos_i(i+3-dim) += fd_step;
		rays[1+i] = std::make_shared<Ray>(
			init_pos_i, init_moment);
	}
}

template <int dim>
void RayBundle<dim>::project_on_analyzer(Vector<3,double> &analyzer) {

	Vector<3,double> &p = rays[0]->cur_moment;

	auto a = (p,p)*analyzer-(analyzer,p)*p;
	a /= a.norm();
	for(int i=0; i<2; i++) {
		cur_transmission_factors[i] *= (a,cur_pols[i]);
		cur_pols[i] = a;
	}
}

template <int dim>
void RayBundle<dim>::simple_geometrical_spreading_update(bool backtrack) {

	// We check that all the rays have the same optical length
	unsigned int i,j;
	for(unsigned int i=0; i<dim; ++i)
		Assert(
			std::abs(rays[i]->cur_optical_length -
				rays[i+1]->cur_optical_length)<1e-6,
			"The rays of a BundleRay do not have the same optical length.");
	
	// We update the cur_jacobian matrix
	for(i=0; i<dim; ++i)
		for(j=0; j<dim; ++j)
			cur_jac(i,j) = 
				( rays[1+j]->cur_pos(i+3-dim)
				- rays[0]->cur_pos(i+3-dim)) / fd_step;

	// We compute the determinant of the cur_jacobian matrix, i.e. the
	// geometrical spreading
	double old_geom_spreading = cur_geom_spreading;
	cur_geom_spreading = compute_geometrical_spreading();

	// We increment the number of crossed caustic if the geometrical
	// spreading has changed its sign or has gone to zero.
	if(old_geom_spreading*cur_geom_spreading<0) {
		std::complex<double> I(0.,1.);
		if(backtrack) {
			cur_transmission_factors[0] *= I;
			cur_transmission_factors[1] *= I;
		}
		else {
			cur_transmission_factors[0] *= -I;
			cur_transmission_factors[0] *= -I;
		}
	}
}

template <int dim>
void RayBundle<dim>::exact_geometrical_spreading_update(
		Vector<3,double> (&pos_derivatives)[dim+1], double step_length) {

	// We check that all the rays have the same optical length
	unsigned int i,j;
	for(unsigned int i=0; i<dim; ++i)
		Assert(
			std::abs(rays[i]->cur_optical_length -
				rays[i+1]->cur_optical_length)<1e-6,
			"The rays of a BundleRay do not have the same optical length.");

	// We save the old value of the geometrical spreading
	for(i=0; i<dim; ++i)
		for(j=0; j<dim; ++j)
			cur_jac(i,j) = 
				( rays[1+j]->cur_pos(i+3-dim)
				- rays[0]->cur_pos(i+3-dim)) / fd_step +
				( pos_derivatives[1+j](i+3-dim) 
				- pos_derivatives[0](i+3-dim)) * (-step_length) / fd_step;
	double old_geom_spreading = compute_geometrical_spreading();
	
	// We compute a slightly shifted version of the gometrical spreading 
	// along the ray
	for(i=0; i<dim; ++i)
		for(j=0; j<dim; ++j)
			cur_jac(i,j) = 
				( rays[1+j]->cur_pos(i+3-dim)
				- rays[0]->cur_pos(i+3-dim)) / fd_step +
				( pos_derivatives[1+j](i+3-dim) 
				- pos_derivatives[0](i+3-dim)) * (fd_step-step_length) / fd_step;
	double shifted_geom_spreading = compute_geometrical_spreading();
	
	// We compute the new value of the geometrical spreading,
	// as well as the derivative of the geometrical spreading with
	// respect to the optical length
	for(i=0; i<dim; ++i)
		for(j=0; j<dim; ++j)
			cur_jac(i,j) = 
				( rays[1+j]->cur_pos(i+3-dim)
				- rays[0]->cur_pos(i+3-dim)) / fd_step;
	cur_geom_spreading = compute_geometrical_spreading();
	double geom_spreading_derivative =
		(shifted_geom_spreading-old_geom_spreading) / fd_step;

	// We fit the variation of the geom spreading with a second-order
	// polynomial, and determine the number of roots between the previous
	// and current update. The origin of the polynomial variable is
	// centered on the old update.
	double a =
		(cur_geom_spreading-old_geom_spreading) / std::pow(step_length,2.) -
		geom_spreading_derivative / step_length;
	double b = geom_spreading_derivative;
	double c = old_geom_spreading;
	double delta = b*b-4*a*c;
	double smin = std::min(0., step_length);
	double smax = std::max(0., step_length);

	int num_roots = 0;
	if(delta==0 && -b/(2*a)>=smin && -b/(2*a)<=smax)
		num_roots++;
	else if(delta>0) {
		if( (-b+std::sqrt(delta))/(2*a)>=smin &&
				(-b+std::sqrt(delta))/(2*a)<=smax)
			num_roots++;
		if( (-b-std::sqrt(delta))/(2*a)>=smin &&
				(-b-std::sqrt(delta))/(2*a)<=smax)
			num_roots++;
	}

	// We increment (or decrement) the number of crossed caustic.
	std::complex<double> I;
	if(step_length<0) {
		cur_transmission_factors[0] *= std::exp(I*num_roots*PI/2.);
		cur_transmission_factors[1] *= std::exp(I*num_roots*PI/2.);
	}
	else {
		cur_transmission_factors[0] *= std::exp(-I*num_roots*PI/2.);
		cur_transmission_factors[1] *= std::exp(-I*num_roots*PI/2.);
	}
}

template <int dim>
double RayBundle<dim>::compute_geometrical_spreading() {

	double geom_spreading;
	if(dim==1)
		geom_spreading =
			cur_jac(0,0);
	else if(dim==2)
		geom_spreading = 
			cur_jac(0,0)*cur_jac(1,1) -
			cur_jac(1,0)*cur_jac(0,1);
	else if(dim==3)
		geom_spreading = 
			cur_jac(0,0)*cur_jac(1,1)*cur_jac(2,2) +
			cur_jac(0,1)*cur_jac(1,2)*cur_jac(2,0) +
			cur_jac(1,0)*cur_jac(0,2)*cur_jac(2,1) -
			cur_jac(2,0)*cur_jac(1,1)*cur_jac(0,2) -
			cur_jac(1,0)*cur_jac(0,1)*cur_jac(2,2) -
			cur_jac(1,2)*cur_jac(2,1)*cur_jac(0,0);
	else
		throw(std::string("Dimension not implemented (RayBundle)."));

	return geom_spreading;
}

template <int dim>
void RayBundle<dim>::compute_E_field(
		Vector<3,std::complex<double> > &res, int pol_idx) {

	Assert(
		pol_idx==0 || pol_idx==1,
		"Error: the polarisation index must be 0 or 1 (RayBundle<dim>::compute_E_field)");

	double a =
		init_rescaled_amplitude /
		( cur_optical_index * std::sqrt(std::abs(cur_geom_spreading)) );
	res = cur_transmission_factors[pol_idx]*a*cur_pols[pol_idx];
}

template <int dim>
void RayBundle<dim>::compute_B_field(
		Vector<3,std::complex<double> > &res, int pol_idx) {

	Assert(
		pol_idx==0 || pol_idx==1,
		"Error: the polarisation index must be 0 or 1 (RayBundle<dim>::compute_B_field)");

	Vector<3,std::complex<double> > E;
	compute_E_field(E, pol_idx);
	res = rays[0]->cur_moment ^ E;
}

template <int dim>
double RayBundle<dim>::optical_length() {
	return rays[0]->cur_optical_length;
}

template <int dim>
bool RayBundle<dim>::is_NA_limited() {

	return NA_limited;
}

template <int dim>
bool RayBundle<dim>::is_evanescent() {

	return evanescent;
}

template class RayBundle<1>;
template class RayBundle<2>;
template class RayBundle<3>;
