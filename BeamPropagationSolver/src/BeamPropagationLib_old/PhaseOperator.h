#ifndef PHASEOPERATOR_H
#define PHASEOPERATOR_H

#include <complex>

#include "PermittivityTensorField.h"

class PhaseOperator {
public:
	PhaseOperator(
		const PermittivityTensorField &eps,
		double wavelength);

	/**
	 * Apply the evolution operator for the phase.
	 */
	void vmult(
		VectorField<std::complex<double> > &dst,
		const VectorField<std::complex<double> > &src) const;

	void z_step_increment();
	void z_step_reset() {
		iz = 0;
	}

private:
	/**
	 * Number of points in each space direction.
	 */
	int Nx, Ny, Nz;

	/**
	 * Wavevector in empty space multiplied by the mesh spacing in the
	 * z-direction
	 */
	const double delta_Z;
	/**
	 * Current z index.
	 */
	int iz;

	/**
	 * A vector field containing the 3 independent components (xx, yy, xy)
	 * of the (half) phase operator exp(I*sqrt(eps_transverse)*delta_Z/2)
	 */
	VectorField<std::complex<double> > op_field;
};

#endif
