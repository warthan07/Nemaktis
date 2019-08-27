#ifndef ADIOPERATORY_H
#define ADIOPERATORY_H

#include "TriDiagonalADIOperator.h"

class ADIOperatorY {
public:
	ADIOperatorY(
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
	class Block00 : public TriDiagonalADIOperatorY {
	public:
		Block00(
			const PermittivityTensorField &eps,
			double wavelength, const BPMSettings& bpm_settings);
	
	protected:
		virtual std::complex<double> diag_inf(int ix, int iy) const;
		virtual std::complex<double> diag_mid(int ix, int iy) const;
		virtual std::complex<double> diag_sup(int ix, int iy) const;
	};
	class Block11 : public TriDiagonalADIOperatorY {
	public:
		Block11(
			const PermittivityTensorField &eps,
			double wavelength, const BPMSettings& bpm_settings);
	
		void vmult_corr(
			VectorField<std::complex<double> > &dst,
			int dst_comp,
			const VectorField<std::complex<double> > &src,
			int src_comp) const;

	protected:
		virtual std::complex<double> diag_inf(int ix, int iy) const;
		virtual std::complex<double> diag_mid(int ix, int iy) const;
		virtual std::complex<double> diag_sup(int ix, int iy) const;
	};
	class Block10 : public BaseADIOperator {
	public:
		Block10(
			const PermittivityTensorField &eps,
			double wavelength, const BPMSettings& bpm_settings);

		virtual void vmult(
			VectorField<std::complex<double> > &dst,
			int dst_comp,
			const VectorField<std::complex<double> > &src,
			int src_comp) const;
	};
	Block00 block_00;
	Block10 block_10;
	Block11 block_11;
	
	int N_woodbury_steps;
};


#endif
