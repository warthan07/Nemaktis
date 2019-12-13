#include "BaseBPMOperator.h"

BaseBPMOperator::BaseBPMOperator(
		const PermittivityTensorField &eps,
		double wavelength, const BPMSettings& bpm_settings) :
	wavevector(2*PI/wavelength),
	delta_X(wavevector*eps.mesh.delta_x),
	delta_Y(wavevector*eps.mesh.delta_y),
	delta_Z(wavevector*eps.mesh.delta_z),
	Nx(eps.mesh.Nx),
	Ny(eps.mesh.Ny),
	Nz(eps.mesh.Nz),
	eps(eps),
	periodic_x(bpm_settings.boundary_types()[0]==BoundaryType::Periodic),
	periodic_y(bpm_settings.boundary_types()[1]==BoundaryType::Periodic),
	iz(0),
	Nz_substeps(1) {}
