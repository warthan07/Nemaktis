#ifndef PHASEEVOLUTIONOPERATOR_H
#define PHASEEVOLUTIONOPERATOR_H

#include "BaseBPMOperator.h"

class PhaseEvolutionOperator : public BaseBPMOperator {
public:
	PhaseEvolutionOperator(
		const PermittivityTensorField& eps,
		double wavelength, const RootSettings &settings);

	/**
	 * Apply this evolution operator to a transverse optical field.
	 */
	virtual int apply(TransverseOpticalField &src);
	/**
	 * Reinit this operator for the current transverse plane.
	 */
	virtual void update();

private:
	/**
	 * Vectors containing the 3 independent components (xx, yy, xy) of the (half) phase
	 * evolution operator exp(mu*sqrt(eps_tr)) with mu=I*delta_Z/2
	 */
	std::vector<std::complex<double> > op_xx, op_yy, op_xy;

	/**
	 * Pure imaginary number I.
	 */
	const std::complex<double> I;

	/**
	 * Reference refractive index, set to the ordinary index of the LC
	 */
	const double no;

    /**
     * Complex factor allowing to take into account a nonzero imaginary part for the
     * extraordinary index, linked to absorption along the director
     */
    const std::complex<double> absorp_fac;	
};

#endif
