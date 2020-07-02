#ifndef PERMITTIVITYTENSORFIELD_H
#define PERMITTIVITYTENSORFIELD_H

#include "VectorField.h"
#include "PhysicsCoefficients.h"

class PermittivityTensorField : public VectorField<double> {
public:
	/**
	 * Initializes the permittivity tensor field from a director or
	 * q-tensor field.
	 */
	PermittivityTensorField(
		const VectorField<double> &lc_sol,
		const PhysicsCoefficients &coefs,
		double wavelength);

	double get_ne() const {
		return ne;
	}
	double get_no() const {
		return no;
	}
	double get_nref() const {
		return nref;
	}

	/**
	 * Get xx component of the permittivity tensor
	 */
	double xx(const Index3D &p) const;
	/**
	 * Get yy component of the permittivity tensor
	 */
	double yy(const Index3D &p) const;
	/**
	 * Get zz component of the permittivity tensor
	 */
	double zz(const Index3D &p) const;
	/**
	 * Get xy component of the permittivity tensor
	 */
	double xy(const Index3D &p) const;
	/**
	 * Get xz component of the permittivity tensor
	 */
	double xz(const Index3D &p) const;
	/**
	 * Get yz component of the permittivity tensor
	 */
	double yz(const Index3D &p) const;
	/**
	 * Get xx component of the transverse permittivity tensor: exx-exz^2/ezz
	 */
	double xx_tr(const Index3D &p) const;
	/**
	 * Get yy component of the transverse permittivity tensor: eyy-eyz^2/ezz
	 */
	double yy_tr(const Index3D &p) const;
	/**
	 * Get xy component of the transverse permittivity tensor: exy-exz*eyz/ezz
	 */
	double xy_tr(const Index3D &p) const;
	/**
	 * Get xx component of the square root of the transverse permittivity tensor
	 */
	double xx_sqrt(const Index3D &p) const;
	/**
	 * Get yy component of the square root of the transverse permittivity tensor
	 */
	double yy_sqrt(const Index3D &p) const;
	/**
	 * Get xy component of the square root of the transverse permittivity tensor
	 */
	double xy_sqrt(const Index3D &p) const;


private:
	/**
	 * Ordinary refractive index sqrt(eps_perp)
	 */
	const double no;
	/**
	 * Extraordinary refractive index sqrt(eps_parallel)
	 */
	const double ne;
	/**
	 * Reference index (mean eigenvalue of the averaged square root transverse permittivity)
	 */
	double nref;
	/**
	 * Relative permittivity in the direction orthogonal to the director.
	 */
	const double eps_perp;
	/**
	 * Anisotropy of relative permittivity.
	 */
	const double eps_a;
	/**
	 * Relative permittivity of the host fluid in case of a non-trivial
	 * domain mask.
	 */
	const double eps_host;
};

#endif
