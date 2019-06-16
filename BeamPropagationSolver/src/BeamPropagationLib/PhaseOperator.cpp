#include <cmath>

#include "PhaseOperator.h"
#include "error.h"

PhaseOperator::PhaseOperator(
		const PermittivityTensorField &eps,
		double wavelength) :
	Nx(eps.mesh.Nx),
	Ny(eps.mesh.Ny),
	Nz(eps.mesh.Nz),
	delta_Z(2*M_PI*eps.mesh.delta_z/wavelength),
	iz(0),
	op_field(eps.mesh, 3) {

	const std::complex<double> I(0,1.);
	const double no = eps.get_no();

	// Parallel loop computing the evolution operator for the phases:
	// exp(I*delta_n*k0_dz/2), where the matrix delta_n is defined for
	// each vertex in the functions delta_n_**.
	#pragma omp parallel for
	for(unsigned int iz=0; iz<Nz-1; iz++) {
		for(unsigned int iy=0; iy<Ny; iy++) {
			for(unsigned int ix=0; ix<Nx; ix++) {
				// Efficient calculation of d=I*(sqrt(eps_transverse)-no*Id)*delta_Z/2
				// with eps_transverse = {{exx-exz^2/ezz,exy-exz*eyz/ezz}
				//                       ,{exy-exz*eyz/ezz,eyy-eyz^2/ezz}}
				// based on its invariants (determinant and half trace)
				double exx = eps.xx({ix,iy,iz});
				double eyy = eps.yy({ix,iy,iz});
				double ezz = eps.zz({ix,iy,iz});
				double exy = eps.xy({ix,iy,iz});
				double exz = eps.xz({ix,iy,iz});
				double eyz = eps.yz({ix,iy,iz});

				double fxx = exx-exz*exz/ezz;
				double fyy = eyy-eyz*eyz/ezz;
				double fxy = exy-exz*eyz/ezz;

				double det = fxx*fyy-fxy*fxy;
				double tr = 0.5*(fxx+fyy);

				// if the eigenvalues are almost equal, tr*tr-det can 
				// be slighly negative due to finite arithmetic error.
				// In this case, we can use a simpler (and numerically
				// stable) formula for the effective index neff
				double neff;
				if(tr*tr-det < 1e-8)
					neff = std::sqrt(tr);
				else
					neff = 
						( std::sqrt(tr+std::sqrt(tr*tr-det)) 
						+ std::sqrt(tr-std::sqrt(tr*tr-det)) ) / 2.;

				double dxx = 0.5*delta_Z*(neff-no+(fxx-fyy)/(4*neff));
				double dyy = 0.5*delta_Z*(neff-no+(fyy-fxx)/(4*neff));
				double dxy = 0.5*delta_Z*fxy/(2*neff);

				
				// Efficient calculation of exp(I*d) based on the
				// Cayleigh-Hamilton theorem and Sylvester formula (see
				// wikipedia page "Matrix exponential").
				double diff = (dxx-dyy)/2.;
				tr = (dxx+dyy)/2.;
				det = std::sqrt(dxy*dxy+diff*diff);
				double sinc_det = (std::abs(det)<1e-8) ? 1 : std::sin(det)/det;

				op_field({ix,iy,iz},0) =
					std::exp(I*tr) * ( std::cos(det) + I*diff*sinc_det );
				op_field({ix,iy,iz},1) =
					std::exp(I*tr) * ( std::cos(det) - I*diff*sinc_det );
				op_field({ix,iy,iz},2) =
					std::exp(I*tr) * I*dxy*sinc_det;
			}
		}
	}
}

void PhaseOperator::vmult(
		VectorField<std::complex<double> > &dst,
		const VectorField<std::complex<double> > &src) const {

	#pragma omp parallel for
	for(unsigned int transverse_idx=0; transverse_idx<Nx*Ny; transverse_idx++) {
		unsigned int ix = (transverse_idx%Nx);
		unsigned int iy = (transverse_idx/Nx);

		dst({ix,iy,0},0) = 
			op_field({ix,iy,iz},0) * src({ix,iy,0},0) +
			op_field({ix,iy,iz},2) * src({ix,iy,0},1);
		dst({ix,iy,0},1) =
			op_field({ix,iy,iz},2) * src({ix,iy,0},0) +
			op_field({ix,iy,iz},1) * src({ix,iy,0},1);
	}
}

void PhaseOperator::z_step_increment() {

	Assert(
		iz<Nz-1, "Out-of-range z increment for the BaseADI operator.");
	iz++;
}
