#ifndef BASEBPMOPERATOR_H
#define BASEBPMOPERATOR_H

#include <set>
#include <complex>
#include <memory>
#include <Eigen/SparseCore>

#include "PhysicsCoefficients.h"
#include "CartesianMesh.h"
#include "PermittivityTensorField.h"
#include "TransparentBoundaryConditions.h"


/**
 * Base class for all BPM operators
 */
class BaseBPMOperator {
public:
	BaseBPMOperator(
		const PermittivityTensorField &eps,
		double wavelength, const BPMSettings& bpm_settings);

	/**
	 * Apply this operator to a transverse optical field.
	 */
	virtual int apply(TransverseOpticalField &src) = 0;
	/**
	 * Reinit this operator for the current transverse plane.
	 */
	virtual void update() = 0;
	/**
	 * Increment the index of the transverse plane
	 */
	void z_step_increment() {
		iz++;
	}
	/**
	 * Reinitialize the index of the transverse plane
	 */
	void z_reinit() {
		iz = 0;
	}
	/**
	 * Increment the number of evolution substeps per z-slab
	 */
	virtual void increment_Nz_substeps() {
		Nz_substeps++;
	}
	/**
	 * Decrement the number of evolution substeps per z-slab
	 */
	virtual void decrement_Nz_substeps() {
		if(Nz_substeps>1)
			Nz_substeps--;
	}

	/**
	 * Wavevector in empty space.
	 */
	const double wavevector;

	/**
	 * Mesh spacings multiplied by the wavevector.
	 */
	const double delta_X, delta_Y, delta_Z;

	/**
	 * Number of points in each space direction.
	 */
	const int Nx, Ny, Nz;

	/**
	 * Reference to the full permittivity tensor field.
	 */
	const PermittivityTensorField& eps;

protected:
	/**
	 * Do we need to apply periodic BC along x?
	 */
	const bool periodic_x;
	/**
	 * Do we need to apply periodic BC along y?
	 */
	const bool periodic_y;

	/**
	 * Current z index.
	 */
	int iz;
	/**
	 * Number of evolution sub-steps per z-slab
	 */
	int Nz_substeps;
};

////////////////////////////////////////////////////////////
// Ugly trick to have a less shitty std::complex on linux //
////////////////////////////////////////////////////////////

#ifndef _MSC_VER

template <typename T>
struct identity_t { typedef T type; };

#define COMPLEX_OPS(OP)                                                                      \
	template <typename _Tp>                                                                  \
	std::complex<_Tp>                                                                        \
	operator OP(std::complex<_Tp> lhs, const typename identity_t<_Tp>::type & rhs) {         \
		return lhs OP rhs;                                                                   \
	}                                                                                        \
	template <typename _Tp>                                                                  \
	std::complex<_Tp>                                                                        \
	operator OP(const typename identity_t<_Tp>::type & lhs, const std::complex<_Tp> & rhs) { \
		return lhs OP rhs;                                                                   \
	}

COMPLEX_OPS(+)
COMPLEX_OPS(-)
COMPLEX_OPS(*)
COMPLEX_OPS(/)

#undef COMPLEX_OPS

#endif

#endif
