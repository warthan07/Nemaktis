#ifndef BASEADIOPERATOR_H
#define BASEADIOPERATOR_H

#include <complex>

#include "PermittivityTensorField.h"
#include "error.h"

class BaseADIOperator {
public:
	BaseADIOperator(
		const PermittivityTensorField &eps,
		double wavelength, const BPMSettings& bpm_settings);

	virtual void vmult(
		VectorField<std::complex<double> > &dst,
		unsigned int dst_comp,
		const VectorField<std::complex<double> > &src,
		unsigned int src_comp) const = 0;

	void switch_to_backward_operator();
	void switch_to_forward_operator();
	void z_step_increment();
	void z_step_reset() {
		iz = 0;
	}

protected:
	std::complex<double> d_ovr_dX(
		const VectorField<std::complex<double> > &src,
		unsigned int ix, unsigned int iy, unsigned int comp) const;
	std::complex<double> d_ovr_dY(
		const VectorField<std::complex<double> > &src,
		unsigned int ix, unsigned int iy, unsigned int comp) const;
	std::complex<double> d2_ovr_dX2(
		const VectorField<std::complex<double> > &src,
		unsigned int ix, unsigned int iy, unsigned int comp) const;
	std::complex<double> d2_ovr_dY2(
		const VectorField<std::complex<double> > &src,
		unsigned int ix, unsigned int iy, unsigned int comp) const;
	std::complex<double> d2_ovr_dXdY(
		const VectorField<std::complex<double> > &src,
		unsigned int ix, unsigned int iy, unsigned int comp) const;

	/**
	 * Wavevector in empty space
	 */
	const double k0;

	/**
	 * Mesh spacings multiplied by the wavevector in empty space.
	 */
	const double delta_X, delta_Y, delta_Z;

	/**
	 * Number of points in each space direction.
	 */
	const unsigned int Nx, Ny, Nz;

	/**
	 * Current z index.
	 */
	unsigned int iz;

	/**
	 * Coefficient used in the expression of the ADI operator. Set to
	 * 1/(2*no)+i*k0*delta_z/2 in the forward mode and
	 * 1/(2*no)-i*k0*delta_z/2 in the backward mode.
	 */
	std::complex<double> mu;
	/**
	 * Imaginary number i.
	 */
	const std::complex<double> I;

	const bool periodic_x;
	const bool periodic_y;

	/**
	 * Reference to the permittivity tensor field object
	 */
	const PermittivityTensorField &eps;
};

///////////////////////////////////////////////////
// Ugly trick to have a less shitty std::complex //
///////////////////////////////////////////////////
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
