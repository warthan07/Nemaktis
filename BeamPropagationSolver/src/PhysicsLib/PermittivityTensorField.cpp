#include "PermittivityTensorField.h"

PermittivityTensorField::PermittivityTensorField(
		const VectorField<double> &lc_sol,
		const PhysicsCoefficients &coefs,
		double wavelength) :
	VectorField<double>(lc_sol.mesh, 12),
	no(coefs.get_no(wavelength)),
	ne(coefs.get_ne(wavelength)),
	eps_perp(std::pow(no,2.)),
	eps_a(std::pow(ne,2.)-eps_perp),
	eps_host(std::pow(coefs.get_nhost(wavelength),2.)) {

	bool mask_val;

	if(lc_sol.field_dim == 3) {
		#pragma omp parallel for firstprivate(mask_val)
		for(int iz=0; iz<mesh.Nz; iz++) {
			for(int iy=0; iy<mesh.Ny; iy++) {
				for(int ix=0; ix<mesh.Nx; ix++) {
					Index3D p{ix,iy,iz};
					mask_val = lc_sol.get_mask_val(p);
					(*this)(p,0) = mask_val ?
						eps_perp + eps_a * std::pow(lc_sol(p,0),2.) : eps_host;
					(*this)(p,1) = mask_val ?
						eps_perp + eps_a * std::pow(lc_sol(p,1),2.) : eps_host;
					(*this)(p,2) = mask_val ?
						eps_perp + eps_a * std::pow(lc_sol(p,2),2.) : eps_host;
					(*this)(p,3) = mask_val ?
						eps_a * lc_sol(p,0)*lc_sol(p,1) : 0;
					(*this)(p,4) = mask_val ?
						eps_a * lc_sol(p,0)*lc_sol(p,2) : 0;
					(*this)(p,5) = mask_val ?
						eps_a * lc_sol(p,1)*lc_sol(p,2) : 0;
				}
			}
		}
	}
	else if(lc_sol.field_dim == 6) {
		double eps_a_eff = 2*eps_a/3;
		double eps_iso = eps_perp+eps_a/3.;
		#pragma omp parallel for firstprivate(mask_val)
		for(int iz=0; iz<mesh.Nz; iz++) {
			for(int iy=0; iy<mesh.Ny; iy++) {
				for(int ix=0; ix<mesh.Nx; ix++) {
					Index3D p{ix,iy,iz};
					mask_val = lc_sol.get_mask_val(p);
					(*this)(p,0) = mask_val ? 
						eps_iso + eps_a_eff * lc_sol(p,0) : eps_host;
					(*this)(p,1) =  mask_val ? 
						eps_iso + eps_a_eff * lc_sol(p,1) : eps_host;
					(*this)(p,2) =  mask_val ? 
						eps_iso + eps_a_eff * lc_sol(p,2) : eps_host;
					(*this)(p,3) =  mask_val ? 
						eps_a_eff * lc_sol(p,3) : 0;
					(*this)(p,4) =  mask_val ? 
						eps_a_eff * lc_sol(p,4) : 0;
					(*this)(p,5) =  mask_val ? 
						eps_a_eff * lc_sol(p,5) : 0;
				}
			}
		}
	}

	// Values of the components of the transverse permittivity tensor and its square root.
	// We also estimate the reference index for the shift operator from the mean of the
	// eigenvalues of the mesh-averaged square root transverse permittivity tensor.
	double det, tr, neff, nref_sum;
	#pragma omp parallel for private(det, tr, neff) reduction(+:nref_sum)
	for(int iz=0; iz<mesh.Nz; iz++) {
		for(int iy=0; iy<mesh.Ny; iy++) {
			for(int ix=0; ix<mesh.Nx; ix++) {
				Index3D p{ix,iy,iz};

				(*this)(p,6) = xx(p) - std::pow(xz(p),2)/zz(p);
				(*this)(p,7) = yy(p) - std::pow(yz(p),2)/zz(p);
				(*this)(p,8) = xy(p) - xz(p)*yz(p)/zz(p);

				det = xx_tr(p)*yy_tr(p)-xy_tr(p)*xy_tr(p);
				tr = 0.5*(xx_tr(p)+yy_tr(p));
				if(tr*tr-det < 1e-8)
					neff = std::sqrt(tr);
				else
					neff = 
						( std::sqrt(tr+std::sqrt(tr*tr-det)) 
						+ std::sqrt(tr-std::sqrt(tr*tr-det)) ) / 2.;
				
				(*this)(p,9)  = neff+(xx_tr(p)-yy_tr(p))/(4*neff);
				(*this)(p,10) = neff+(yy_tr(p)-xx_tr(p))/(4*neff);
				(*this)(p,11) = xy_tr(p)/(2*neff);

				nref_sum += neff;
			}
		}
	}
	nref = nref_sum/(mesh.Nx*mesh.Ny*mesh.Nz);
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
double PermittivityTensorField::xx_tr(const Index3D &p) const {
	return (*this)(p, 6);
}
double PermittivityTensorField::yy_tr(const Index3D &p) const {
	return (*this)(p, 7);
}
double PermittivityTensorField::xy_tr(const Index3D &p) const {
	return (*this)(p, 8);
}
double PermittivityTensorField::xx_sqrt(const Index3D &p) const {
	return (*this)(p, 9);
}
double PermittivityTensorField::yy_sqrt(const Index3D &p) const {
	return (*this)(p, 10);
}
double PermittivityTensorField::xy_sqrt(const Index3D &p) const {
	return (*this)(p, 11);
}
