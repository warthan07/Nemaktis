#ifndef TRIDIAGONALADIOPERATOR_H
#define TRIDIAGONALADIOPERATOR_H

#include "BaseADIOperator.h"

class TriDiagonalADIOperator : public BaseADIOperator {
public:
	TriDiagonalADIOperator(
		const PermittivityTensorField &eps,
		double wavelength, const BPMSettings& bpm_settings);

protected:
	virtual std::complex<double> diag_inf(unsigned int ix, unsigned int iy) const = 0;
	virtual std::complex<double> diag_mid(unsigned int ix, unsigned int iy) const = 0;
	virtual std::complex<double> diag_sup(unsigned int ix, unsigned int iy) const = 0;
};

class TriDiagonalADIOperatorX : public TriDiagonalADIOperator {
public:
	TriDiagonalADIOperatorX(
		const PermittivityTensorField &eps,
		double wavelength, const BPMSettings& bpm_settings);

	virtual void vmult(
		VectorField<std::complex<double> > &dst,
		unsigned int dst_comp,
		const VectorField<std::complex<double> > &src,
		unsigned int src_comp) const;
	void inverse_vmult(
		VectorField<std::complex<double> > &dst,
		unsigned int dst_comp,
		const VectorField<std::complex<double> > &src,
		unsigned int src_comp) const;

	void update_tbc_wavectors(
		const VectorField<std::complex<double> > &src,
		unsigned int comp);

protected:
	/**
	 * Array storing the transparent boundary condition factors
	 * exp(i*kx*delta_x) for the lower X boundary
	 */
	std::vector<std::complex<double> > upper_tbc_exp;
	/**
	 * Array storing the transparent boundary condition factors
	 * exp(i*kx*delta_x) for the upper X boundary
	 */
	std::vector<std::complex<double> > lower_tbc_exp;
};

class TriDiagonalADIOperatorY : public TriDiagonalADIOperator {
public:
	TriDiagonalADIOperatorY(
		const PermittivityTensorField &eps,
		double wavelength, const BPMSettings& bpm_settings);

	virtual void vmult(
		VectorField<std::complex<double> > &dst,
		unsigned int dst_comp,
		const VectorField<std::complex<double> > &src,
		unsigned int src_comp) const;
	void inverse_vmult(
		VectorField<std::complex<double> > &dst,
		unsigned int dst_comp,
		const VectorField<std::complex<double> > &src,
		unsigned int src_comp) const;

	void update_tbc_wavectors(
		const VectorField<std::complex<double> > &src,
		unsigned int comp);

protected:
	/**
	 * Array storing the transparent boundary condition factors
	 * exp(i*ky*delta_y) for the left Y boundary
	 */
	std::vector<std::complex<double> > left_tbc_exp;
	/**
	 * Array storing the transparent boundary condition factors
	 * exp(i*ky*delta_y) for the right Y boundary
	 */
	std::vector<std::complex<double> > right_tbc_exp;
};

#endif
