#ifndef FRESNELOPERATOR_H
#define FRESNELOPERATOR_H

#include <complex>
#include "PermittivityTensorField.h"

class FresnelOperator {
public: 
	FresnelOperator(
		const PermittivityTensorField &eps,
		double input_refractive_index);

	/**
	 * Apply the fresnel boundary conditions at the input interface.
	 */
	void vmult(
		VectorField<std::complex<double> > &dst,
		const VectorField<std::complex<double> > &src) const;
private:
	/**
	 * Number of points in each space direction.
	 */
	int Nx, Ny;
	/**
	 * Input refractive index.
	 */
	const double ni;
	
	/**
	 * Reference to the permittivity tensor field
	 */
	const PermittivityTensorField &eps;
};

#endif
