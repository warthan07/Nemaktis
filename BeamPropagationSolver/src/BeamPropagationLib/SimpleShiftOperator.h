#ifndef SIMPLESHIFTOPERATOR_H
#define SIMPLESHIFTOPERATOR_H

#include "BaseBPMOperator.h"
#include "BlockSparseMatrix.h"

class SimpleShiftOperator : public BaseBPMOperator {
public:
	SimpleShiftOperator(
		const PermittivityTensorField& eps, double wavelength,
		std::pair<double,double> q_val, const RootSettings &settings);

	/**
	 * Apply this operator to a transverse optical field.
	 */
	virtual int apply(TransverseOpticalField &src);
	/**
	 * This operator does not need to be updated since it only depends on nref.
	 */
	virtual void update() {};

private:
	inline int xwrap(int ix) {
		if(ix>Nx-2)
			return ix%(Nx-1);
		else if(ix<0)
			return ix+(1-ix/(Nx-1))*(Nx-1);
		else
			return ix;
	}
	inline int ywrap(int iy) {
		if(iy>Ny-2)
			return iy%(Ny-1);
		else if(iy<0)
			return iy+(1-iy/(Ny-1))*(Ny-1);
		else
			return iy;
	}

	/**
	 * Polynomial weights used for cubic interpolation
	 */
	double pol_weight_x[4], pol_weight_y[4];
	/**
	 * Shifted index for the translation operators
	 */
	std::vector<int> shifted_ix[4], shifted_iy[4];
	/**
	 * Unitary complex number applying a phase correction
	 */
	std::complex<double> phase_factor;

	/**
	 * A temporary transverse field used in the apply method.
	 */
	TransverseOpticalField tmp;
};


#endif
