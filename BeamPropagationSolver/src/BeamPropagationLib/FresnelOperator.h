#ifndef FRESNELOPERATOR_H
#define FRESNELOPERATOR_H

#include <complex>
#include "PermittivityTensorField.h"
#include "OpticalField.h"

class FresnelOperator {
public: 
	FresnelOperator(
		const PermittivityTensorField &eps,
		double refractive_index);

	/**
	 * Apply the fresnel boundary conditions at the interface.
	 */
	virtual void apply(TransverseOpticalField &src) const = 0;

protected:
	/**
	 * Number of points in each space direction.
	 */
	int Nx, Ny;
	/**
	 * Refractive index of isotropic medium.
	 */
	const double niso;
	
	/**
	 * Reference to the permittivity tensor field
	 */
	const PermittivityTensorField &eps;
};


class InputFresnelOperator : public FresnelOperator {
public: 
	InputFresnelOperator(
		const PermittivityTensorField &eps,
		double input_refractive_index);

	/**
	 * Apply the fresnel boundary conditions at the input interface.
	 */
	virtual void apply(TransverseOpticalField &src) const;
};

class OutputFresnelOperator : public FresnelOperator {
public: 
	OutputFresnelOperator(
		const PermittivityTensorField &eps,
		double output_refractive_index);

	/**
	 * Apply the fresnel boundary conditions at the input interface.
	 */
	virtual void apply(TransverseOpticalField &src) const;
};


#endif
