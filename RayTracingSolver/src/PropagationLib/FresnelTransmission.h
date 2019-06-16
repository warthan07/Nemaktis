#ifndef FRESNELTRANSMISSION_H
#define FRESNELTRANSMISSION_H

#include <iostream>
#include <Eigen/Dense>

#include "ODEFunction.h"

template <int dim>
class FresnelTransmission {
public:
	FresnelTransmission() {};

	void update(
		std::shared_ptr<ODEFunction<dim> > (&ode_functions)[4],
		RayBundle<dim> &ray_bundle);
	void update(
		std::shared_ptr<ODEFunction<dim> > (&ode_functions)[4],
		RayBundle<dim> &ray_bundle_1,
		RayBundle<dim> &ray_bundle_2);

private:
	void compute_fresnel_coefs(
		std::shared_ptr<ODEFunction<dim> > (&ode_functions)[4],
		Ray &incident_ray, Vector<3,std::complex<double> > &incident_pol);

	Eigen::Matrix4cd fresnel_matrix;
	Eigen::Vector4cd fresnel_rhs, fresnel_coefs;

	Vector<3,std::complex<double> > e_pol[4];
	Vector<3,std::complex<double> > b_pol[4];
	bool is_evanescent[4];

	Vector<3,std::complex<double> > e_i;
	Vector<3,std::complex<double> > b_i;

	Vector<3,double> t_s;
	Vector<3,double> t_p;
};

#endif
