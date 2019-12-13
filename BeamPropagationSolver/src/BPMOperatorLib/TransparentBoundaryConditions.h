#ifndef TRANSPARENTBOUNDARYCONDITIONS_H
#define TRANSPARENTBOUNDARYCONDITIONS_H

#include <complex>

#include "OpticalField.h"
#include "PhysicsCoefficients.h"

class TransparentBoundaryConditions {
public:
	TransparentBoundaryConditions(int Nx, int Ny);

	/**
	 * Update the TBC coefs from the values of the transverse field for the current transverse
	 * plane.
	 */
	void update_coefs(const TransverseOpticalField &transverse_field);

	/**
	 * Distribute the TBC constraints on an Eigen vector.
	 */
	void apply_constraints(Eigen::VectorXcd &src) const;

	/**
	 * Get a TBC coef for the +X or -X boundary (no OOB check, be careful!)
	 */
	std::complex<double> coef_plusX(int iy, int field_comp) const {
		return _tbc_coef_plusX[field_comp][iy-1];
	}
	std::complex<double> coef_minusX(int iy, int field_comp) const {
		return _tbc_coef_minusX[field_comp][iy-1];
	}

	/**
	 * Get a TBC coef for the +Y or -Y boundary (no OOB check, be careful!)
	 */
	std::complex<double> coef_plusY(int ix, int field_comp) const {
		return _tbc_coef_plusY[field_comp][ix-1];
	}
	std::complex<double> coef_minusY(int ix, int field_comp) const {
		return _tbc_coef_minusY[field_comp][ix-1];
	}

	/**
	 * Number of points in the x and y space direction, and total number
	 * of points in a transverse plane N=Nx*Ny.
	 */
	const int Nx, Ny, N;

private:
	/**
	 * TBC coefs for the +X and -X boundaries
	 */
	std::vector<std::complex<double> > _tbc_coef_plusX[2];
	std::vector<std::complex<double> > _tbc_coef_minusX[2];

	/**
	 * TBC coefs for the +Y and -Y boundaries
	 */
	std::vector<std::complex<double> > _tbc_coef_plusY[2];
	std::vector<std::complex<double> > _tbc_coef_minusY[2];
};

#endif
