#include "FresnelOperator.h"

FresnelOperator::FresnelOperator(
		const PermittivityTensorField &eps,
		double input_refractive_index) :
	Nx(eps.mesh.Nx),
	Ny(eps.mesh.Ny),
	eps(eps),
	ni(input_refractive_index) {}

void FresnelOperator::vmult(
		VectorField<std::complex<double> > &dst,
		const VectorField<std::complex<double> > &src) const {

	const std::complex<double> I(0,1.);

	// Parallel loop which applies the Fresnel boudary conditions to the input fields.
	#pragma omp parallel for
	for(unsigned int iperp=0; iperp<Nx*Ny; iperp++) {
		unsigned int ix = iperp%Nx;
		unsigned int iy = iperp/Nx;

		// Efficient calculation of d=sqrt(eps_transverse)
		// with eps_transverse = {{exx-exz^2/ezz,exy-exz*eyz/ezz}
		//                       ,{exy-exz*eyz/ezz,eyy-eyz^2/ezz}}
		// based on its invariants (determinant and half trace)
		double exx = eps.xx({ix,iy,0});
		double eyy = eps.yy({ix,iy,0});
		double ezz = eps.zz({ix,iy,0});
		double exy = eps.xy({ix,iy,0});
		double exz = eps.xz({ix,iy,0});
		double eyz = eps.yz({ix,iy,0});

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

		double dxx = (neff+(fxx-fyy)/(4*neff));
		double dyy = (neff+(fyy-fxx)/(4*neff));
		double dxy = fxy/(2*neff);


		// We apply the Fresnel BC
		det = (ni+dxx)*(ni+dyy)-dxy*dxy;
		dst({ix,iy,0},0) = 
			2*ni*(ni+dyy)/det * src({ix,iy,0},0) -
			2*ni*dxy/det * src({ix,iy,0},1);
		dst({ix,iy,0},1) =
			 - 2*ni*dxy/det * src({ix,iy,0},0) +
			 2*ni*(ni+dxx)/det * src({ix,iy,0},1);
	}
}
