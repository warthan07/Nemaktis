#include "FresnelOperator.h"

FresnelOperator::FresnelOperator(
		const PermittivityTensorField &eps,
		double input_refractive_index) :
	Nx(eps.mesh.Nx),
	Ny(eps.mesh.Ny),
	eps(eps),
	ni(input_refractive_index) {}

void FresnelOperator::apply(TransverseOpticalField &src) const {

	const std::complex<double> I(0,1.);

	// Parallel loop which applies the Fresnel boudary conditions to the input fields.
	#pragma omp parallel for
	for(int iperp=0; iperp<Nx*Ny; iperp++) {
		int ix = iperp%Nx;
		int iy = iperp/Nx;

		double nxx = eps.xx_sqrt({ix,iy,0});
		double nyy = eps.yy_sqrt({ix,iy,0});
		double nxy = eps.xy_sqrt({ix,iy,0});

		// We apply the Fresnel BC
		double det = (ni+nxx)*(ni+nyy)-nxy*nxy;
		std::complex<double> new_field_val[2] = {
			2*ni*(ni+nyy)/det * src({ix,iy},0) - 2*ni*nxy/det * src({ix,iy},1),
			- 2*ni*nxy/det * src({ix,iy},0) + 2*ni*(ni+nxx)/det * src({ix,iy},1)};

		src({ix,iy},0) = new_field_val[0];
		src({ix,iy},1) = new_field_val[1];
	}
}
