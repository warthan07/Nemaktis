#include "PhaseEvolutionOperator.h"

PhaseEvolutionOperator::PhaseEvolutionOperator(
		const PermittivityTensorField& eps,
		double wavelength, const RootSettings &settings) :
	BaseBPMOperator(eps, wavelength),
	op_xx(Nx*Ny),
	op_yy(Nx*Ny),
	op_xy(Nx*Ny),
	I(0,1),
	no(eps.get_no()),
    absorp_fac(std::complex<double>(1, eps.get_ne_imag()/(eps.get_ne()-no))) {}

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

	#pragma omp parallel for
	for(int iy=0; iy<Ny; iy++) {
		for(int ix=0; ix<Nx; ix++) {
			Index3D p{ix,iy,iz};
			Index3D p_pz{ix,iy,iz+1};

			double w = eps.get_interp_weight(p, 0);
			double w_pz = eps.get_interp_weight(p, 1);

			// Components of tensor d = (sqrt(eps_tr)-no)*delta_Z/2
			double dxx = 0.5*delta_Z*(w*eps.xx_sqrt(p)+w_pz*eps.xx_sqrt(p_pz)-no);
			double dyy = 0.5*delta_Z*(w*eps.yy_sqrt(p)+w_pz*eps.yy_sqrt(p_pz)-no);
			double dxy = 0.5*delta_Z*(w*eps.xy_sqrt(p)+w_pz*eps.xy_sqrt(p_pz));

			// Efficient calculation of exp(I*absorp_fac*d) based on the
			// Cayleigh-Hamilton theorem and Sylvester formula (see
			// wikipedia page "Matrix exponential").
			double diff = (dxx-dyy)/2.;
			double tr = (dxx+dyy)/2.;
			double det = std::sqrt(dxy*dxy+diff*diff);
            std::complex<double> sinc_det = (std::abs(det)<1e-8) ? 1 : std::sin(absorp_fac*det)/det;

			op_xx[ix+Nx*iy] = std::exp(I*absorp_fac*tr) * ( std::cos(absorp_fac*det) + I*diff*sinc_det );
			op_yy[ix+Nx*iy] = std::exp(I*absorp_fac*tr) * ( std::cos(absorp_fac*det) - I*diff*sinc_det );
			op_xy[ix+Nx*iy] = std::exp(I*absorp_fac*tr) * I*dxy*sinc_det;
		}
	}
}
