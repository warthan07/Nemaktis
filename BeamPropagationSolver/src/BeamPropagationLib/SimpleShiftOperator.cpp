#include "SimpleShiftOperator.h"

SimpleShiftOperator::SimpleShiftOperator(
		const PermittivityTensorField& eps, double wavelength,
		std::pair<double,double> q_val, const RootSettings &settings) :
	BaseBPMOperator(eps, wavelength),
	tmp(eps.mesh) {

	double qx = q_val.first;
	double qy = q_val.second;
	double nref = eps.get_nref();

	phase_factor = std::exp(-std::complex<double>(0,delta_Z*(qx*qx+qy*qy)/(4*nref)));

	if(Nx>3) {
		double xshift = -delta_Z*qx/(2*nref*delta_X); // global shift
		double mesh_xshift = std::floor(xshift); // integer shift
		double nu = xshift-mesh_xshift; // local noninteger shift inside a cell

		for(int ix=0; ix<Nx; ix++) {
			shifted_ix[0].push_back(xwrap(ix+mesh_xshift-1));
			shifted_ix[1].push_back(xwrap(ix+mesh_xshift));
			shifted_ix[2].push_back(xwrap(ix+mesh_xshift+1));
			shifted_ix[3].push_back(xwrap(ix+mesh_xshift+2));
		}

		pol_weight_x[0] = -nu*std::pow(1-nu,2)/2;
		pol_weight_x[1] = (1-nu)*(1+nu-3*nu*nu/2);
		pol_weight_x[2] = nu*(1+4*nu-3*nu*nu)/2;
		pol_weight_x[3] = -nu*nu*(1-nu)/2;
	}

	if(Ny>3) {
		double yshift = -delta_Z*qy/(2*nref*delta_X); // global shift
		double mesh_yshift = std::floor(yshift); // integer shift
		double nu = yshift-mesh_yshift; // local noninteger shift inside a cell

		for(int iy=0; iy<Ny; iy++) {
			shifted_iy[0].push_back(ywrap(iy+mesh_yshift-1));
			shifted_iy[1].push_back(ywrap(iy+mesh_yshift));
			shifted_iy[2].push_back(ywrap(iy+mesh_yshift+1));
			shifted_iy[3].push_back(ywrap(iy+mesh_yshift+2));
		}

		pol_weight_y[0] = -nu*std::pow(1-nu,2)/2;
		pol_weight_y[1] = (1-nu)*(1+nu-3*nu*nu/2);
		pol_weight_y[2] = nu*(1+4*nu-3*nu*nu)/2;
		pol_weight_y[3] = -nu*nu*(1-nu)/2;
	}
}

int SimpleShiftOperator::apply(TransverseOpticalField &src) {

	if(Nx>3) {
		#pragma omp parallel for
		for(int iy=0; iy<Ny; iy++) {
			for(int comp=0; comp<2; comp++) {
				for(int ix=0; ix<Nx; ix++) {
					tmp({ix,iy},comp) = 0;
					for(int k=0; k<4; k++)
						tmp({ix,iy},comp) +=
							pol_weight_x[k]*src({shifted_ix[k][ix],iy},comp);
				}
			}
		}
	}
	else
		tmp = src;

	if(Ny>3) {
		#pragma omp parallel for
		for(int ix=0; ix<Nx; ix++) {
			for(int comp=0; comp<2; comp++) {
				for(int iy=0; iy<Ny; iy++) {
					src({ix,iy},comp) = 0;
					for(int k=0; k<4; k++)
						src({ix,iy},comp) +=
							pol_weight_y[k]*tmp({ix,shifted_iy[k][iy]},comp);
				}
			}
		}
	}
	else
		src = tmp;

	#pragma omp parallel for
	for(int ix=0; ix<Nx; ix++)
		for(int comp=0; comp<2; comp++)
			for(int iy=0; iy<Ny; iy++)
				src({ix,iy},comp) *= phase_factor;

	return 1;
}
