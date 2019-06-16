#include "TriDiagonalADIOperator.h"

TriDiagonalADIOperator::TriDiagonalADIOperator(
		const PermittivityTensorField &eps,
		double wavelength, const BPMSettings& bpm_settings) :
	BaseADIOperator(eps, wavelength, bpm_settings) {}

TriDiagonalADIOperatorX::TriDiagonalADIOperatorX(
		const PermittivityTensorField &eps,
		double wavelength, const BPMSettings& bpm_settings) :
	TriDiagonalADIOperator(eps, wavelength, bpm_settings) {

	if(bpm_settings.wide_angle_corrections()[0])
		mu.real(1./(2*eps.get_no()));

	if(periodic_x) {
		upper_tbc_exp.assign(Ny,0);
		lower_tbc_exp.assign(Ny,0);
	}
	else {
		upper_tbc_exp.assign(Ny,1);
		lower_tbc_exp.assign(Ny,1);
	}
}


void TriDiagonalADIOperatorX::vmult(
		VectorField<std::complex<double> > &dst, 
		unsigned int dst_comp,
		const VectorField<std::complex<double> > &src,
		unsigned int src_comp) const {

	if(periodic_x) {
		#pragma omp parallel for
		for(unsigned int iy=0; iy<Ny; iy++) {
			dst({0,iy,0},dst_comp) =
				diag_inf(0,iy) * src({Nx-1,iy,0},src_comp) +
				diag_mid(0,iy) * src({0,iy,0},src_comp) +
				diag_sup(0,iy) * src({1,iy,0},src_comp);
			dst({Nx-1,iy,0},dst_comp) =
				diag_inf(Nx-1,iy) * src({Nx-2,iy,0},src_comp) +
				diag_mid(Nx-1,iy) * src({Nx-1,iy,0},src_comp) +
				diag_sup(Nx-1,iy) * src({0,iy,0},src_comp);
			for(unsigned int ix=1; ix<Nx-1; ix++) {
				dst({ix,iy,0},dst_comp) =
					diag_inf(ix,iy) * src({ix-1,iy,0},src_comp) +
					diag_mid(ix,iy) * src({ix,iy,0},src_comp) +
					diag_sup(ix,iy) * src({ix+1,iy,0},src_comp);
			}
		}
	}
	else {
		#pragma omp parallel for
		for(unsigned int iy=0; iy<Ny; iy++) {
			dst({0,iy,0},dst_comp) =
				diag_mid(0,iy) * src({0,iy,0},src_comp) +
				diag_sup(0,iy) * src({1,iy,0},src_comp);
			dst({Nx-1,iy,0},dst_comp) =
				diag_inf(Nx-1,iy) * src({Nx-2,iy,0},src_comp) +
				diag_mid(Nx-1,iy) * src({Nx-1,iy,0},src_comp);
			for(unsigned int ix=1; ix<Nx-1; ix++) {
				dst({ix,iy,0},dst_comp) =
					diag_inf(ix,iy) * src({ix-1,iy,0},src_comp) +
					diag_mid(ix,iy) * src({ix,iy,0},src_comp) +
					diag_sup(ix,iy) * src({ix+1,iy,0},src_comp);
			}
		}
	}
}

void TriDiagonalADIOperatorX::inverse_vmult(
		VectorField<std::complex<double> > &dst, 
		unsigned int dst_comp,
		const VectorField<std::complex<double> > &src,
		unsigned int src_comp) const {

	if(periodic_x) {
		#pragma omp parallel for
		for(unsigned int iy=0; iy<Ny; iy++) {
			std::vector<std::complex<double> > gamma(Nx-1);
			std::vector<std::complex<double> > corr(Nx);

			std::complex<double> beta = 2.*diag_mid(0,iy);
			dst({0,iy,0},dst_comp) = src({0,iy,0},src_comp) / beta;
			corr[0] = -0.5;
			
			for(unsigned int ix=1; ix<Nx-1; ix++) {
				gamma[ix-1] = diag_sup(ix-1,iy) / beta;
				beta = diag_mid(ix,iy) - diag_inf(ix,iy)*gamma[ix-1];

				dst({ix,iy,0},dst_comp) =
					( src({ix,iy,0},src_comp) - diag_inf(ix,iy)*dst({ix-1,iy,0},dst_comp) ) / beta;
				corr[ix] = - diag_inf(ix,iy)*corr[ix-1] / beta;
			}

			gamma[Nx-2] = diag_sup(Nx-2,iy) / beta;
			beta =
				diag_mid(Nx-1,iy) + diag_sup(Nx-1,iy)*diag_inf(0,iy)/diag_mid(0,iy)
				- diag_inf(Nx-1,iy)*gamma[Nx-2];

			dst({Nx-1,iy,0},dst_comp) =
				( src({Nx-1,iy,0},src_comp) - diag_inf(Nx-1,iy)*dst({Nx-2,iy,0},dst_comp) ) / beta;
			corr[Nx-1] = ( diag_sup(Nx-1,iy) - diag_inf(Nx-1,iy)*corr[Nx-2] ) / beta;

			for(unsigned int ix=Nx-1; ix>0; ix--) {
				dst({ix-1,iy,0},dst_comp) -= gamma[ix-1]*dst({ix,iy,0},dst_comp);
				corr[ix-1] -= gamma[ix-1]*corr[ix];
			}

			std::complex<double> frac =
				( dst({0,iy,0},dst_comp) - diag_inf(0,iy)/diag_mid(0,iy)*dst({Nx-1,iy,0},dst_comp) ) /
				( 1 + corr[0] - diag_inf(0,iy)/diag_mid(0,iy)*corr[Nx-1] );
			for(unsigned int ix=0; ix<Nx; ix++)
				dst({ix,iy,0},dst_comp) -= frac*corr[ix];
		}
	}
	else {
		#pragma omp parallel for
		for(unsigned int iy=0; iy<Ny; iy++) {
			std::vector<std::complex<double> > gamma(Nx-1);
			std::complex<double> beta = diag_mid(0, iy);
			dst({0,iy,0},dst_comp) = src({0,iy,0},src_comp) / beta;
			
			for(unsigned int ix=1; ix<Nx; ix++) {
				gamma[ix-1] = diag_sup(ix-1, iy) / beta;
				beta = diag_mid(ix, iy) - diag_inf(ix,iy)*gamma[ix-1];
				dst({ix,iy,0},dst_comp) =
					( src({ix,iy,0},src_comp) - diag_inf(ix,iy)*dst({ix-1,iy,0},dst_comp) ) / beta;
			}
			for(unsigned int ix=Nx-1; ix>0; ix--)
				dst({ix-1,iy,0},dst_comp) -= gamma[ix-1]*dst({ix,iy,0},dst_comp);
		}
	}
}

void TriDiagonalADIOperatorX::update_tbc_wavectors(
		const VectorField<std::complex<double> > &src,
		unsigned int comp) {

	if(!periodic_x) {
		#pragma omp parallel for
		for(unsigned int iy=0; iy<Ny; iy++) {
			if(std::abs(src({1,iy,0},comp))!=0)
				lower_tbc_exp[iy] = src({0,iy,0},comp)/src({1,iy,0},comp);
			else
				lower_tbc_exp[iy] = 1;
			if(std::imag(lower_tbc_exp[iy])<0)
				lower_tbc_exp[iy] = std::real(lower_tbc_exp[iy]);
	
			if(std::abs(src({Nx-2,iy,0},comp))!=0)
				upper_tbc_exp[iy] = src({Nx-1,iy,0},comp)/src({Nx-2,iy,0},comp);
			else
				upper_tbc_exp[iy] = 1;
			if(std::imag(upper_tbc_exp[iy])<0)
				upper_tbc_exp[iy] = std::real(lower_tbc_exp[iy]);
		}
	}
}


TriDiagonalADIOperatorY::TriDiagonalADIOperatorY(
		const PermittivityTensorField &eps,
		double wavelength, const BPMSettings& bpm_settings) :
	TriDiagonalADIOperator(eps, wavelength, bpm_settings) {

	if(bpm_settings.wide_angle_corrections()[1])
		mu.real(1./(2*eps.get_no()));

	if(periodic_y) {
		left_tbc_exp.assign(Nx, 0.),
		right_tbc_exp.assign(Nx, 0.);
	}
	else {
		left_tbc_exp.assign(Nx, 1.),
		right_tbc_exp.assign(Nx, 1.);
	}
}

void TriDiagonalADIOperatorY::vmult(
		VectorField<std::complex<double> > &dst,
		unsigned int dst_comp,
		const VectorField<std::complex<double> > &src,
		unsigned int src_comp) const {

	if(periodic_y) {
		#pragma omp parallel for
		for(unsigned int ix=0; ix<Nx; ix++) {
			dst({ix,0,0},dst_comp) =
				diag_inf(ix,0) * src({ix,Ny-1,0},src_comp) +
				diag_mid(ix,0) * src({ix,0,0},src_comp) +
				diag_sup(ix,0) * src({ix,1,0},src_comp);
			dst({ix,Ny-1,0},dst_comp) =
				diag_inf(ix,Ny-1) * src({ix,Ny-2,0},src_comp) +
				diag_mid(ix,Ny-1) * src({ix,Ny-1,0},src_comp) +
				diag_sup(ix,Ny-1) * src({ix,0,0},src_comp);
			for(unsigned int iy=1; iy<Ny-1; iy++) {
				dst({ix,iy,0},dst_comp) =
					diag_inf(ix, iy) * src({ix,iy-1,0},src_comp) +
					diag_mid(ix, iy) * src({ix,iy,0},src_comp) +
					diag_sup(ix, iy) * src({ix,iy+1,0},src_comp);
			}
		}
	}
	else {
		#pragma omp parallel for
		for(unsigned int ix=0; ix<Nx; ix++) {
			dst({ix,0,0},dst_comp) =
				diag_mid(ix,0) * src({ix,0,0},src_comp) +
				diag_sup(ix,0) * src({ix,1,0},src_comp);
			dst({ix,Ny-1,0},dst_comp) =
				diag_inf(ix,Ny-1) * src({ix,Ny-2,0},src_comp) +
				diag_mid(ix,Ny-1) * src({ix,Ny-1,0},src_comp);
			for(unsigned int iy=1; iy<Ny-1; iy++) {
				dst({ix,iy,0},dst_comp) =
					diag_inf(ix, iy) * src({ix,iy-1,0},src_comp) +
					diag_mid(ix, iy) * src({ix,iy,0},src_comp) +
					diag_sup(ix, iy) * src({ix,iy+1,0},src_comp);
			}
		}
	}
}

void TriDiagonalADIOperatorY::inverse_vmult(
		VectorField<std::complex<double> > &dst,
		unsigned int dst_comp,
		const VectorField<std::complex<double> > &src,
		unsigned int src_comp) const {

	if(periodic_y) {
		#pragma omp parallel for
		for(unsigned int ix=0; ix<Nx; ix++) {
			std::vector<std::complex<double> > gamma(Ny-1);
			std::vector<std::complex<double> > corr(Ny);

			std::complex<double> beta = 2.*diag_mid(ix,0);
			dst({ix,0,0},dst_comp) = src({ix,0,0},src_comp) / beta;
			corr[0] = -0.5;
			
			for(unsigned int iy=1; iy<Ny-1; iy++) {
				gamma[iy-1] = diag_sup(ix,iy-1) / beta;
				beta = diag_mid(ix,iy) - diag_inf(ix,iy)*gamma[iy-1];

				dst({ix,iy,0},dst_comp) =
					( src({ix,iy,0},src_comp) - diag_inf(ix,iy)*dst({ix,iy-1,0},dst_comp)) / beta;
				corr[iy] = - diag_inf(ix,iy)*corr[iy-1] / beta;
			}

			gamma[Ny-2] = diag_sup(ix,Ny-2) / beta;
			beta =
				diag_mid(ix,Ny-1) + diag_sup(ix,Ny-1)*diag_inf(ix,0)/diag_mid(ix,0)
				- diag_inf(ix,Ny-1)*gamma[Ny-2];

			dst({ix,Ny-1,0},dst_comp) =
				( src({ix,Ny-1,0},src_comp) - diag_inf(ix,Ny-1)*dst({ix,Ny-2,0},dst_comp) ) / beta;
			corr[Ny-1] = ( diag_sup(ix,Ny-1) - diag_inf(ix,Ny-1)*corr[Ny-2] ) / beta;

			for(unsigned int iy=Ny-1; iy>0; iy--) {
				dst({ix,iy-1,0},dst_comp) -= gamma[iy-1]*dst({ix,iy,0},dst_comp);
				corr[iy-1] -= gamma[iy-1]*corr[iy];
			}

			std::complex<double> frac =
				( dst({ix,0,0},dst_comp) - diag_inf(ix,0)/diag_mid(ix,0)*dst({ix,Ny-1,0},dst_comp) ) /
				( 1 + corr[0] - diag_inf(ix,0)/diag_mid(ix,0)*corr[Ny-1] );
			for(unsigned int iy=0; iy<Ny; iy++)
				dst({ix,iy,0},dst_comp) -= frac*corr[iy];
		}
	}
	else {
		#pragma omp parallel for
		for(unsigned int ix=0; ix<Nx; ix++) {
			std::vector<std::complex<double> > gamma(Ny-1);
			std::complex<double> beta = diag_mid(ix, 0);
			dst({ix,0,0},dst_comp) = src({ix,0,0},src_comp) / beta;
			
			for(unsigned int iy=1; iy<Ny; iy++) {
				gamma[iy-1] = diag_sup(ix,iy-1) / beta;
				beta = diag_mid(ix,iy) - diag_inf(ix,iy)*gamma[iy-1];
				dst({ix,iy,0},dst_comp) =
					( src({ix,iy,0},src_comp) - diag_inf(ix,iy)*dst({ix,iy-1,0},dst_comp) ) / beta;
			}
			for(unsigned int iy=Ny-1; iy>0; iy--)
				dst({ix,iy-1,0},dst_comp) -= gamma[iy-1]*dst({ix,iy,0},dst_comp);
		}
	}
}

void TriDiagonalADIOperatorY::update_tbc_wavectors(
		const VectorField<std::complex<double> > &src,
		unsigned int comp) {

	if(!periodic_y) {
		#pragma omp parallel for
		for(unsigned int ix=0; ix<Nx; ix++) {
			if(std::abs(src({ix,1,0},comp))!=0)
				left_tbc_exp[ix] = src({ix,0,0},comp)/src({ix,1,0},comp);
			else
				left_tbc_exp[ix] = 1;
			if(std::imag(left_tbc_exp[ix])<0)
				left_tbc_exp[ix] = std::real(left_tbc_exp[ix]);
	
			if(std::abs(src({ix,Ny-2,0},comp))!=0)
				right_tbc_exp[ix] = src({ix,Ny-1,0},comp)/src({ix,Ny-2,0},comp);
			else
				right_tbc_exp[ix] = 1;
			if(std::imag(right_tbc_exp[ix])<0)
				right_tbc_exp[ix] = std::real(left_tbc_exp[ix]);
		}
	}
}
