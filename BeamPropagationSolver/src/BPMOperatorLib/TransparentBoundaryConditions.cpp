#include "TransparentBoundaryConditions.h"

TransparentBoundaryConditions::TransparentBoundaryConditions(
		int Nx, int Ny) :
	Nx(Nx),
	Ny(Ny),
	N(Nx*Ny) {
	
	for(int comp=0; comp<2; comp++) {
		_tbc_coef_plusX[comp].resize(Ny-2,1.);
		_tbc_coef_minusX[comp].resize(Ny-2,1.);
		_tbc_coef_plusY[comp].resize(Nx-2,1.);
		_tbc_coef_minusY[comp].resize(Nx-2,1.);
	}
}

void TransparentBoundaryConditions::update_coefs(
		const TransverseOpticalField &src) {

	#pragma omp parallel for
	for(int iy=1; iy<Ny-1; iy++) {
		for(int comp=0; comp<2; comp++) {
			if(std::abs(src({2,iy},comp))!=0)
				_tbc_coef_minusX[comp][iy-1] = src({1,iy},comp)/src({2,iy},comp);
			else
				_tbc_coef_minusX[comp][iy-1] = 1;
			if(std::imag(_tbc_coef_minusX[comp][iy-1])<0)
				_tbc_coef_minusX[comp][iy-1] = std::real(_tbc_coef_minusX[comp][iy-1]);

			if(std::abs(src({Nx-3,iy},comp))!=0)
				_tbc_coef_plusX[comp][iy-1] = src({Nx-2,iy},comp)/src({Nx-3,iy},comp);
			else
				_tbc_coef_plusX[comp][iy-1] = 1;
			if(std::imag(_tbc_coef_plusX[comp][iy-1])<0)
				_tbc_coef_plusX[comp][iy-1] = std::real(_tbc_coef_plusX[comp][iy-1]);
		}
	}

	#pragma omp parallel for
	for(int ix=1; ix<Nx-1; ix++) {
		for(int comp=0; comp<2; comp++) {
			if(std::abs(src({ix,2},comp))!=0)
				_tbc_coef_minusY[comp][ix-1] = src({ix,1},comp)/src({ix,2},comp);
			else
				_tbc_coef_minusY[comp][ix-1] = 1;
			if(std::imag(_tbc_coef_minusY[comp][ix-1])<0)
				_tbc_coef_minusY[comp][ix-1] = std::real(_tbc_coef_minusY[comp][ix-1]);

			if(std::abs(src({ix,Ny-3},comp))!=0)
				_tbc_coef_plusY[comp][ix-1] = src({ix,Ny-2},comp)/src({ix,Ny-3},comp);
			else
				_tbc_coef_plusY[comp][ix-1] = 1;
			if(std::imag(_tbc_coef_plusY[comp][ix-1])<0)
				_tbc_coef_plusY[comp][ix-1] = std::real(_tbc_coef_plusY[comp][ix-1]);
		}
	}
}

void TransparentBoundaryConditions::apply_constraints(Eigen::VectorXcd &src) const {

	#pragma omp parallel for
	for(int iy=1; iy<Ny-1; iy++) {
		for(int comp=0; comp<2; comp++) {
			src[Nx*iy+comp*N] = src[1+Nx*iy+comp*N]*coef_minusX(iy,comp);
			src[Nx-1+Nx*iy+comp*N] = src[Nx-2+Nx*iy+comp*N]*coef_plusX(iy,comp);
		}
	}

	#pragma omp parallel for
	for(int ix=1; ix<Nx-1; ix++) {
		for(int comp=0; comp<2; comp++) {
			src[ix+comp*N] = src[ix+Nx+comp*N]*coef_minusY(ix,comp);
			src[ix+Nx*(Ny-1)+comp*N] = src[ix+Nx*(Ny-2)+comp*N]*coef_plusY(ix,comp);
		}
	}

	for(int comp=0; comp<2; comp++) {
		src[comp*N] = src[Nx+1+comp*N]*coef_minusX(1,comp)*coef_minusY(1,comp);
		src[Nx-1+comp*N] = src[2*Nx-2+comp*N]*coef_plusX(1,comp)*coef_minusY(Nx-2,comp);
		src[Nx*(Ny-1)+comp*N] = src[1+Nx*(Ny-2)+comp*N]*coef_minusX(Ny-2,comp)*coef_plusY(1,comp);
		src[N-1+comp*N] = src[N-Nx-2+comp*N]*coef_plusX(Ny-2,comp)*coef_plusY(Nx-2,comp);
	}
}
