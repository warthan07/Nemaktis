#include "PermittivityTensorField_old.h"

PermittivityTensorField::PermittivityTensorField(
		const VectorField<double> &lc_sol,
		const PhysicsCoefficients &coefs,
		double wavelength) :
	VectorField<double>(lc_sol.mesh, 6),
	no(coefs.get_no(wavelength)),
	ne(coefs.get_ne(wavelength)),
	eps_perp(std::pow(no,2.)),
	eps_a(std::pow(ne,2.)-eps_perp),
	eps_host(std::pow(coefs.get_nhost(wavelength),2.)) {

	bool mask_val0, mask_val1;

	if(lc_sol.field_dim == 3) {
		#pragma omp parallel for firstprivate(mask_val0, mask_val1)
		for(int iz=0; iz<mesh.Nz-1; iz++) {
			for(int iy=0; iy<mesh.Ny; iy++) {
				for(int ix=0; ix<mesh.Nx; ix++) {
					mask_val0 = lc_sol.get_mask_val({ix,iy,iz});
					mask_val1 = lc_sol.get_mask_val({ix,iy,iz+1});
	
					(*this)({ix,iy,iz},0) = 0.5 * (
						(mask_val0 ?
						 	eps_perp+eps_a*std::pow(lc_sol({ix,iy,iz},0),2.) : 	eps_host) +
						(mask_val1 ?
						 	eps_perp+eps_a*std::pow(lc_sol({ix,iy,iz+1},0),2.) : eps_host) );
					(*this)({ix,iy,iz},1) = 0.5 * (
						(mask_val0 ?
						 	eps_perp+eps_a*std::pow(lc_sol({ix,iy,iz},1),2.) : eps_host) +
						(mask_val1 ?
						 	eps_perp+eps_a*std::pow(lc_sol({ix,iy,iz+1},1),2.) : eps_host) );
					(*this)({ix,iy,iz},2) = 0.5 * (
						(mask_val0 ?
						 	eps_perp+eps_a*std::pow(lc_sol({ix,iy,iz},2),2.) : eps_host) +
						(mask_val1 ?
						 	eps_perp+eps_a*std::pow(lc_sol({ix,iy,iz+1},2),2.) : eps_host) );
					(*this)({ix,iy,iz},3) = 0.5 * (
						(mask_val0 ?
						 	eps_a*lc_sol({ix,iy,iz},0)*lc_sol({ix,iy,iz},1) : 0) +
						(mask_val1 ?
						 	eps_a*lc_sol({ix,iy,iz+1},0)*lc_sol({ix,iy,iz+1},1) : 0) );
					(*this)({ix,iy,iz},4) = 0.5 * (
						(mask_val0 ?
						 	eps_a*lc_sol({ix,iy,iz},0)*lc_sol({ix,iy,iz},2) : 0) +
						(mask_val1 ?
						 	eps_a*lc_sol({ix,iy,iz+1},0)*lc_sol({ix,iy,iz+1},2) : 0) );
					(*this)({ix,iy,iz},5) = 0.5 * (
						(mask_val0 ?
						 	eps_a*lc_sol({ix,iy,iz},1)*lc_sol({ix,iy,iz},2) : 0) +
						(mask_val1 ?
						 	eps_a*lc_sol({ix,iy,iz+1},1)*lc_sol({ix,iy,iz+1},2) : 0) );
				}
			}
		}
	}
	else if(lc_sol.field_dim == 6) {
		double eps_a_eff = 2*eps_a/3;
		double eps_iso = eps_perp+eps_a/3.;
		#pragma omp parallel for firstprivate(mask_val0, mask_val1)
		for(int iz=0; iz<mesh.Nz-1; iz++) {
			for(int iy=0; iy<mesh.Ny; iy++) {
				for(int ix=0; ix<mesh.Nx; ix++) {
					mask_val0 = lc_sol.get_mask_val({ix,iy,iz});
					mask_val1 = lc_sol.get_mask_val({ix,iy,iz+1});
	
					(*this)({ix,iy,iz},0) = 0.5 * (
						(mask_val0 ? eps_iso+eps_a_eff*lc_sol({ix,iy,iz},0) : eps_host) +
						(mask_val1 ? eps_iso+eps_a_eff*lc_sol({ix,iy,iz+1},0) : eps_host) );
					(*this)({ix,iy,iz},1) = 0.5 * (
						(mask_val0 ? eps_iso+eps_a_eff*lc_sol({ix,iy,iz},1) : eps_host) +
						(mask_val1 ? eps_iso+eps_a_eff*lc_sol({ix,iy,iz+1},1) : eps_host) );
					(*this)({ix,iy,iz},2) = 0.5 * (
						(mask_val0 ? eps_iso+eps_a_eff*lc_sol({ix,iy,iz},2) : eps_host) +
						(mask_val1 ? eps_iso+eps_a_eff*lc_sol({ix,iy,iz+1},2) : eps_host) );
					(*this)({ix,iy,iz},3) = 0.5 * (
						(mask_val0 ? eps_a_eff*lc_sol({ix,iy,iz},3) : 0) +
						(mask_val1 ? eps_a_eff*lc_sol({ix,iy,iz+1},3) : 0) );
					(*this)({ix,iy,iz},4) = 0.5 * (
						(mask_val0 ? eps_a_eff*lc_sol({ix,iy,iz},4) : 0) +
						(mask_val1 ? eps_a_eff*lc_sol({ix,iy,iz+1},4) : 0) );
					(*this)({ix,iy,iz},5) = 0.5 * (
						(mask_val0 ? eps_a_eff*lc_sol({ix,iy,iz},5) : 0) +
						(mask_val1 ? eps_a_eff*lc_sol({ix,iy,iz+1},5) : 0) );
				}
			}
		}
	}
}

double PermittivityTensorField::xx(const Index3D &p) const {
	return (*this)(p, 0);
}
double PermittivityTensorField::yy(const Index3D &p) const {
	return (*this)(p, 1);
}
double PermittivityTensorField::zz(const Index3D &p) const {
	return (*this)(p, 2);
}
double PermittivityTensorField::xy(const Index3D &p) const {
	return (*this)(p, 3);
}
double PermittivityTensorField::xz(const Index3D &p) const {
	return (*this)(p, 4);
}
double PermittivityTensorField::yz(const Index3D &p) const {
	return (*this)(p, 5);
}
