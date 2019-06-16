#ifndef ADIOPERATORX_H
#define ADIOPERATORX_H

#include "TriDiagonalADIOperator.h"

class ADIOperatorX {
public:
	ADIOperatorX(
		const PermittivityTensorField &eps,
		double wavelength, const BPMSettings& bpm_settings);

	void vmult(
		VectorField<std::complex<double> > &dst,
		const VectorField<std::complex<double> > &src) const;
	void inverse_vmult(
		VectorField<std::complex<double> > &dst,
		const VectorField<std::complex<double> > &src) const;

	void switch_to_backward_operator();
	void switch_to_forward_operator();
	void z_step_increment();
	void z_step_reset();
	void update_tbc_wavectors(const VectorField<std::complex<double> > &src);

private:
	class Block00 : public TriDiagonalADIOperatorX {
	public:
		Block00(
			const PermittivityTensorField &eps,
			double wavelength, const BPMSettings& bpm_settings);

		void vmult_corr(
			VectorField<std::complex<double> > &dst,
			unsigned int dst_comp,
			const VectorField<std::complex<double> > &src,
			unsigned int src_comp) const;
	
	protected:
		virtual std::complex<double> diag_inf(unsigned int ix, unsigned int iy) const;
		virtual std::complex<double> diag_mid(unsigned int ix, unsigned int iy) const;
		virtual std::complex<double> diag_sup(unsigned int ix, unsigned int iy) const;
	};
	class Block11 : public TriDiagonalADIOperatorX {
	public:
		Block11(
			const PermittivityTensorField &eps,
			double wavelength, const BPMSettings& bpm_settings);
	
	protected:
		virtual std::complex<double> diag_inf(unsigned int ix, unsigned int iy) const;
		virtual std::complex<double> diag_mid(unsigned int ix, unsigned int iy) const;
		virtual std::complex<double> diag_sup(unsigned int ix, unsigned int iy) const;
	};
	class Block01 : public BaseADIOperator {
	public:
		Block01(
			const PermittivityTensorField &eps,
			double wavelength, const BPMSettings& bpm_settings);

		virtual void vmult(
			VectorField<std::complex<double> > &dst,
			unsigned int dst_comp,
			const VectorField<std::complex<double> > &src,
			unsigned int src_comp) const;
	};

	Block00 block_00;
	Block01 block_01;
	Block11 block_11;

	const unsigned int N_woodbury_steps;
};

#endif
