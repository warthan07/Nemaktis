#include "BaseADIOperator.h"

#define M_PI 3.1415926535897932

BaseADIOperator::BaseADIOperator(
		const PermittivityTensorField &eps,
		double wavelength, const BPMSettings& bpm_settings) :
	k0(2*M_PI/wavelength),
	delta_X(k0*eps.mesh.delta_x),
	delta_Y(k0*eps.mesh.delta_y),
	delta_Z(k0*eps.mesh.delta_z),
	Nx(eps.mesh.Nx),
	Ny(eps.mesh.Ny),
	Nz(eps.mesh.Nz),
	iz(0),
	mu(std::complex<double>(0, delta_Z/2)),
	I(std::complex<double>(0,1)),
	periodic_x(bpm_settings.boundary_types()[0]==BoundaryType::Periodic),
	periodic_y(bpm_settings.boundary_types()[1]==BoundaryType::Periodic),
	eps(eps) {}

void BaseADIOperator::switch_to_backward_operator() {
	if(mu.imag()>0)
		mu = std::conj(mu);
}
void BaseADIOperator::switch_to_forward_operator() {
	if(mu.imag()<0)
		mu = std::conj(mu);
}
void BaseADIOperator::z_step_increment() {

	iz++;
}

std::complex<double> BaseADIOperator::d_ovr_dX(
		const VectorField<std::complex<double> > &src,
		int ix, int iy, int comp) const {

	int ix1, ix2;

	if(periodic_x) {
		ix1 = (ix!=0) ? ix-1 : Nx-2;
		ix2 = (ix!=Nx-1) ? ix+1 : 1;
		return
			( src({ix2,iy,0},comp) - src({ix1,iy,0},comp) ) / (2*delta_X);
	}
	else {
		if(ix==0) 			{ ix1 = ix;		ix2 = ix+1;	}
		else if(ix==Nx-1) 	{ ix1 = ix-1;	ix2 = ix;	}
		else				{ ix1 = ix-1;	ix2 = ix+1;	}
	
		return 
			( src({ix2,iy,0},comp) - src({ix1,iy,0},comp) ) / ((ix2-ix1)*delta_X);
	}
}

std::complex<double> BaseADIOperator::d_ovr_dY(
		const VectorField<std::complex<double> > &src,
		int ix, int iy, int comp) const {

	int iy1, iy2;

	if(periodic_y) {
		iy1 = (iy!=0) ? iy-1 : Ny-2; 	
		iy2 = (iy!=Ny-1) ? iy+1 : 1;

		return
			( src({ix,iy2,0},comp) - src({ix,iy1,0},comp) ) / (2*delta_Y);
	}
	else {
		if(iy==0)			{ iy1 = iy;		iy2 = iy+1; }
		else if(iy==Ny-1)	{ iy1 = iy-1;	iy2 = iy;	}
		else 				{ iy1 = iy-1;	iy2 = iy+1; }
	
		return 
			( src({ix,iy2,0},comp) - src({ix,iy1,0},comp) ) / ((iy2-iy1)*delta_Y);
	}
}

std::complex<double> BaseADIOperator::d2_ovr_dX2(
		const VectorField<std::complex<double> > &src,
		int ix, int iy, int comp) const {

	int ix1, ix2, ix3;

	if(periodic_x) {
		ix1 = (ix!=0) ? ix-1 : Nx-2;
		ix3 = (ix!=Nx-1) ? ix+1 : 1;

		return 
			( src({ix1,iy,0},comp) - 2.*src({ix,iy,0},comp) + src({ix3,iy,0},comp) )
			/ (delta_X*delta_X);
	}
	else {
		if(ix==0) 			{	ix1 = ix;		ix2 = ix+1;		ix3 = ix+2;	}
		else if(ix==Nx-1) 	{	ix1 = ix-2;		ix2 = ix-1;		ix3 = ix;	}
		else				{ 	ix1 = ix-1;		ix2 = ix;		ix3 = ix+1;	}
	
		return 
			( src({ix1,iy,0},comp) - 2.*src({ix2,iy,0},comp) + src({ix3,iy,0},comp) )
			/ (delta_X*delta_X);
	}
}

std::complex<double> BaseADIOperator::d2_ovr_dY2(
		const VectorField<std::complex<double> > &src,
		int ix, int iy, int comp) const {

	int iy1, iy2, iy3;

	if(periodic_y) {
		iy1 = (iy!=0) ? iy-1 : Ny-2; 	
		iy3 = (iy!=Ny-1) ? iy+1 : 1;

		return 
			( src({ix,iy1,0},comp) - 2.*src({ix,iy,0},comp) + src({ix,iy3,0},comp) )
			/ (delta_Y*delta_Y);
	}
	else {
		if(iy==0) 			{	iy1 = iy;		iy2 = iy+1;		iy3 = iy+2;	}
		else if(iy==Ny-1) 	{	iy1 = iy-2;		iy2 = iy-1;		iy3 = iy;	}
		else				{ 	iy1 = iy-1;		iy2 = iy;		iy3 = iy+1;	}
	
		return 
			( src({ix,iy1,0},comp) - 2.*src({ix,iy2,0},comp) + src({ix,iy3,0},comp) )
			/ (delta_Y*delta_Y);
	}
}

std::complex<double> BaseADIOperator::d2_ovr_dXdY(
		const VectorField<std::complex<double> > &src,
		int ix, int iy, int comp) const {

	int ix1, ix2, iy1, iy2;
	double dX, dY;
	if(periodic_x) {
		ix1 = (ix!=0) ? ix-1 : Nx-2;
		ix2 = (ix!=Nx-1) ? ix+1 : 1;
		dX = 2*delta_X;
	}
	else {
		if(ix==0) 			{ ix1 = ix;		ix2 = ix+1;	}
		else if(ix==Nx-1) 	{ ix1 = ix-1;	ix2 = ix;	}
		else				{ ix1 = ix-1;	ix2 = ix+1;	}
		dX = (ix2-ix1)*delta_X;
	}

	if(periodic_x) {
		iy1 = (iy!=0) ? iy-1 : Ny-2; 	
		iy2 = (iy!=Ny-1) ? iy+1 : 1;
		dY = 2*delta_Y;
	}
	else {
		if(iy==0)			{ iy1 = iy;		iy2 = iy+1; }
		else if(iy==Ny-1)	{ iy1 = iy-1;	iy2 = iy;	}
		else 				{ iy1 = iy-1;	iy2 = iy+1; }
		dY = (iy2-iy1)*delta_Y;
	}

	return 
		( src({ix1,iy1,0},comp) - src({ix1,iy2,0},comp)
		+ src({ix2,iy2,0},comp) - src({ix2,iy1,0},comp) ) / ( dX*dY );
}
