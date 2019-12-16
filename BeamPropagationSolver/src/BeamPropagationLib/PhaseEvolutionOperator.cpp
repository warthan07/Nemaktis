#include "PhaseEvolutionOperator.h"

PhaseEvolutionOperator::PhaseEvolutionOperator(
		const PermittivityTensorField& eps,
		double wavelength, const RootSettings &settings) :
	BaseBPMOperator(eps, wavelength, settings.algorithm.bpm),
	op_xx(Nx*Ny),
	op_yy(Nx*Ny),
	op_xy(Nx*Ny),
	I(0,1),
	nref(eps.get_no()) {}

int PhaseEvolutionOperator::apply(TransverseOpticalField &field) {

	std::complex<double> val_x, val_y; 

	#pragma omp parallel for private(val_x,val_y)
	for(int mesh_idx=0; mesh_idx<Nx*Ny; mesh_idx++) {
		val_x =
			op_xx[mesh_idx] * field(mesh_idx,0) +
			op_xy[mesh_idx] * field(mesh_idx,1);
		val_y =
			op_xy[mesh_idx] * field(mesh_idx,0) +
			op_yy[mesh_idx] * field(mesh_idx,1);

		field(mesh_idx,0) = val_x;
		field(mesh_idx,1) = val_y;
	}

	return 1;
}

void PhaseEvolutionOperator::update() {

	// Parallel loop computing the phase evolution operator inside the LC layer:
	// exp(mu*(sqrt(eps_tr)-nref)).
	#pragma omp parallel for
	for(int iy=0; iy<Ny; iy++) {
		for(int ix=0; ix<Nx; ix++) {
			Index3D p({ix,iy,iz});
			Index3D p_pz({ix,iy,iz+1});

			// Components of tensor d = (sqrt(eps_tr)-nref)*delta_Z/2
			double dxx = 0.25*delta_Z*(eps.xx_sqrt(p)+eps.xx_sqrt(p_pz)-2*nref);
			double dyy = 0.25*delta_Z*(eps.yy_sqrt(p)+eps.yy_sqrt(p_pz)-2*nref);
			double dxy = 0.25*delta_Z*(eps.xy_sqrt(p)+eps.xy_sqrt(p_pz));
			
			// Efficient calculation of exp(I*d) based on the
			// Cayleigh-Hamilton theorem and Sylvester formula (see
			// wikipedia page "Matrix exponential").
			double diff = (dxx-dyy)/2.;
			double tr = (dxx+dyy)/2.;
			double det = std::sqrt(dxy*dxy+diff*diff);
			double sinc_det = (std::abs(det)<1e-8) ? 1 : std::sin(det)/det;

			op_xx[ix+Nx*iy] = std::exp(I*tr) * ( std::cos(det) + I*diff*sinc_det );
			op_yy[ix+Nx*iy] = std::exp(I*tr) * ( std::cos(det) - I*diff*sinc_det );
			op_xy[ix+Nx*iy] = std::exp(I*tr) * I*dxy*sinc_det;
		}
	}
	
}
