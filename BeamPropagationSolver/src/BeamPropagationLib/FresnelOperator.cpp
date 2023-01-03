#include "FresnelOperator.h"

FresnelOperator::FresnelOperator(
		const PermittivityTensorField &eps,
		double refractive_index) :
	Nx(eps.mesh.Nx),
	Ny(eps.mesh.Ny),
	eps(eps),
	niso(refractive_index) {}

InputFresnelOperator::InputFresnelOperator(
		const PermittivityTensorField &eps,
		double input_refractive_index) :
	FresnelOperator(eps, input_refractive_index) {}

void InputFresnelOperator::apply(TransverseOpticalField &src) const {

	const std::complex<double> I(0,1.);

	// Parallel loop which applies the Fresnel boudary conditions
	// at the input interface
	#pragma omp parallel for
	for(int iperp=0; iperp<Nx*Ny; iperp++) {
		int ix = iperp%Nx;
		int iy = iperp/Nx;

		double nxx = eps.xx_sqrt({ix,iy,0});
		double nyy = eps.yy_sqrt({ix,iy,0});
		double nxy = eps.xy_sqrt({ix,iy,0});

		// We apply the Fresnel BC
		double det = (niso+nxx)*(niso+nyy)-nxy*nxy;
		std::complex<double> new_field_val[2] = {
			2*niso*(niso+nyy)/det * src({ix,iy},0) - 2*niso*nxy/det * src({ix,iy},1),
			- 2*niso*nxy/det * src({ix,iy},0) + 2*niso*(niso+nxx)/det * src({ix,iy},1)};

		src({ix,iy},0) = new_field_val[0];
		src({ix,iy},1) = new_field_val[1];
	}
}

OutputFresnelOperator::OutputFresnelOperator(
		const PermittivityTensorField &eps,
		double output_refractive_index) :
	FresnelOperator(eps, output_refractive_index) {}

void OutputFresnelOperator::apply(TransverseOpticalField &src) const {

	const std::complex<double> I(0,1.);
	int Nz = eps.mesh.Nz;

	// Parallel loop which applies the Fresnel boudary conditions
	// at the output interface
	#pragma omp parallel for
	for(int iperp=0; iperp<Nx*Ny; iperp++) {
		int ix = iperp%Nx;
		int iy = iperp/Nx;

		double nxx = eps.xx_sqrt({ix,iy,Nz-1});
		double nyy = eps.yy_sqrt({ix,iy,Nz-1});
		double nxy = eps.xy_sqrt({ix,iy,Nz-1});

		// We apply the Fresnel BC
		double det = (niso+nxx)*(niso+nyy)-nxy*nxy;
		std::complex<double> new_field_val[2] = {
			2.*( (1.-niso*(niso+nyy)/det) * src({ix,iy},0) + niso*nxy/det * src({ix,iy},1) ),
			2.*( (1.-niso*(niso+nxx)/det) * src({ix,iy},1) + niso*nxy/det * src({ix,iy},0) )};

		src({ix,iy},0) = new_field_val[0];
		src({ix,iy},1) = new_field_val[1];
	}
}
