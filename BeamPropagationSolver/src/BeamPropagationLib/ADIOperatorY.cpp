#include "ADIOperatorY.h"

ADIOperatorY::ADIOperatorY(
		const PermittivityTensorField &eps,
		double wavelength, const BPMSettings& bpm_settings) :
	block_00(eps, wavelength, bpm_settings),
	block_10(eps, wavelength, bpm_settings),
	block_11(eps, wavelength, bpm_settings),
	N_woodbury_steps(bpm_settings.n_woodbury_steps) {}

void ADIOperatorY::vmult(
		VectorField<std::complex<double> > &dst,
		const VectorField<std::complex<double> > &src) const {

    VectorField<std::complex<double> > tmp_1(src.mesh, 1, 0);
    VectorField<std::complex<double> > tmp_2(src.mesh, 1, 0);

    block_00.vmult(dst, 0, src, 0);
    block_11.vmult(dst, 1, src, 1);
    block_10.vmult(tmp_1, 0, src, 0);
    block_11.vmult_corr(tmp_2, 0, src, 1);

    #pragma omp parallel for
    for(unsigned int iy=0; iy<src.mesh.Ny; iy++)
        for(unsigned int ix=0; ix<src.mesh.Nx; ix++)
            dst({ix,iy,0},1) += tmp_1({ix,iy,0},0) + tmp_2({ix,iy,0},0);
}

void ADIOperatorY::inverse_vmult(
		VectorField<std::complex<double> > &dst,
		const VectorField<std::complex<double> > &src) const {

	VectorField<std::complex<double> > tmp_1(src.mesh, 1, 0);
	VectorField<std::complex<double> > tmp_2(src.mesh, 1, 0);

	block_00.inverse_vmult(dst, 0, src, 0);
	block_10.vmult(tmp_1, 0, dst, 0);

	#pragma omp parallel for
	for(unsigned int iy=0; iy<src.mesh.Ny; iy++)
		for(unsigned int ix=0; ix<src.mesh.Nx; ix++)
			tmp_1({ix,iy,0},0) = src({ix,iy,0},1)-tmp_1({ix,iy,0},0);

	block_11.inverse_vmult(dst, 1, tmp_1, 0);

	// We use the Woodbury formula to take into account the correction
	// block
	#pragma omp parallel for
	for(unsigned int iy=0; iy<src.mesh.Ny; iy++)
		for(unsigned int ix=0; ix<src.mesh.Nx; ix++)
			tmp_1({ix,iy,0},0) = dst({ix,iy,0},1);
	for(unsigned int it=0; it<N_woodbury_steps; it++) {
		block_11.vmult_corr(tmp_2, 0, tmp_1, 0);
		block_11.inverse_vmult(tmp_1, 0, tmp_2, 0);
		#pragma omp parallel for
		for(unsigned int iy=0; iy<src.mesh.Ny; iy++)
			for(unsigned int ix=0; ix<src.mesh.Nx; ix++)
				dst({ix,iy,0},1) += tmp_1({ix,iy,0},0);
	}
}

void ADIOperatorY::switch_to_backward_operator() {

	block_00.switch_to_backward_operator();
	block_11.switch_to_backward_operator();
	block_10.switch_to_backward_operator();
}

void ADIOperatorY::switch_to_forward_operator() {

	block_00.switch_to_forward_operator();
	block_11.switch_to_forward_operator();
	block_10.switch_to_forward_operator();
}

void ADIOperatorY::z_step_increment() {

	block_00.z_step_increment();
	block_11.z_step_increment();
	block_10.z_step_increment();
}

void ADIOperatorY::z_step_reset() {

	block_00.z_step_reset();
	block_11.z_step_reset();
	block_10.z_step_reset();
}

void ADIOperatorY::update_tbc_wavectors(
		const VectorField<std::complex<double> > &src) {

	block_00.update_tbc_wavectors(src, 0);
	block_11.update_tbc_wavectors(src, 1);
}

ADIOperatorY::Block00::Block00(
		const PermittivityTensorField &eps,
		double wavelength, const BPMSettings& bpm_settings) :
	TriDiagonalADIOperatorY(eps, wavelength, bpm_settings) { 

	if(bpm_settings.wide_angle_corrections()[1])
		mu.real(1./(2*eps.get_no()));
}

std::complex<double> ADIOperatorY::Block00::diag_inf(
		unsigned int ix, unsigned int iy) const {

	double exx = eps.xx({ix,iy,iz});
	return
		mu / (2.*std::sqrt(exx)*std::pow(delta_Y, 2.));
}

std::complex<double> ADIOperatorY::Block00::diag_mid(
		unsigned int ix, unsigned int iy) const {

	double exx = eps.xx({ix,iy,iz});
	if(iy!=0 && iy!=Ny-1)
		return
			1 - mu / (std::sqrt(exx)*std::pow(delta_Y, 2.));
	else if(iy==0)
		return
			1 + mu * (left_tbc_exp[ix]-2) / (2.*std::sqrt(exx)*std::pow(delta_Y, 2.));
	else if(iy==Ny-1)
		return
			1 + mu * (right_tbc_exp[ix]-2) / (2.*std::sqrt(exx)*std::pow(delta_Y, 2.));
	else
		return 0;
}

std::complex<double> ADIOperatorY::Block00::diag_sup(
		unsigned int ix, unsigned int iy) const {

	double exx = eps.xx({ix,iy,iz});
	return
		mu / (2.*std::sqrt(exx)*std::pow(delta_Y, 2.));
}

ADIOperatorY::Block11::Block11(
		const PermittivityTensorField &eps,
		double wavelength, const BPMSettings& bpm_settings) :
	TriDiagonalADIOperatorY(eps, wavelength, bpm_settings) { 

	if(bpm_settings.wide_angle_corrections()[1])
		mu.real(1./(2*eps.get_no()));
}

std::complex<double> ADIOperatorY::Block11::diag_inf(
		unsigned int ix, unsigned int iy) const {

	double eyy = eps.yy({ix,iy,iz});
	double eyz = eps.yz({ix,iy,iz});
	double ezz = eps.zz({ix,iy,iz});
	return
		1./(2*ezz) * (
			mu*std::sqrt(eyy)/std::pow(delta_Y, 2.) +
			I*std::conj(mu)*eyz/delta_Y );
}

std::complex<double> ADIOperatorY::Block11::diag_mid(
		unsigned int ix, unsigned int iy) const {

	double eyy = eps.yy({ix,iy,iz});
	double ezz = eps.zz({ix,iy,iz});
	if(iy!=0 && iy!=Ny-1)
		return
			1 - mu*std::sqrt(eyy)/(ezz*std::pow(delta_Y, 2.));
	else if(iy==0)
		return
			1 + mu*std::sqrt(eyy)*(left_tbc_exp[iy]-2)/(2*ezz*std::pow(delta_Y,2.));
	else if(iy==Ny-1)
		return
			1 + mu*std::sqrt(eyy)*(right_tbc_exp[iy]-2)/(2*ezz*std::pow(delta_Y,2.));
	else
		return 0;
}

std::complex<double> ADIOperatorY::Block11::diag_sup(
		unsigned int ix, unsigned int iy) const {

	double eyy = eps.yy({ix,iy,iz});
	double eyz = eps.yz({ix,iy,iz});
	double ezz = eps.zz({ix,iy,iz});
	return
		1./(2*ezz) * (
			mu*std::sqrt(eyy)/std::pow(delta_Y, 2.) -
			I*std::conj(mu)*eyz/delta_Y );
}

void ADIOperatorY::Block11::vmult_corr(
		VectorField<std::complex<double> > &dst,
		unsigned int dst_comp,
		const VectorField<std::complex<double> > &src,
		unsigned int src_comp) const {

	double eyy, ezz, exy;
	#pragma omp parallel for firstprivate(eyy,ezz,exy)
	for(unsigned int iy=0; iy<Ny; iy++) {
		for(unsigned int ix=0; ix<Nx; ix++) {
			eyy = eps.yy({ix,iy,iz});
			ezz = eps.zz({ix,iy,iz});
			exy = eps.xy({ix,iy,iz});
			dst({ix,iy,0},dst_comp) = 
				mu*exy/(2*ezz*std::sqrt(eyy)) *
				d2_ovr_dXdY(src, ix, iy, src_comp);
		}
	}
}

ADIOperatorY::Block10::Block10(
		const PermittivityTensorField &eps,
		double wavelength, const BPMSettings& bpm_settings) :
	BaseADIOperator(eps, wavelength, bpm_settings) { 

	if(bpm_settings.wide_angle_corrections()[1])
		mu.real(1./(2*eps.get_no()));
}

void ADIOperatorY::Block10::vmult(
		VectorField<std::complex<double> > &dst,
		unsigned int dst_comp,
		const VectorField<std::complex<double> > &src,
		unsigned int src_comp) const {

	double exx, eyy, ezz, exy, exz, eyz;
	#pragma omp parallel for firstprivate(exx,eyy,ezz,exy,exz,eyz)
	for(unsigned int iy=0; iy<Ny; iy++) {
		for(unsigned int ix=0; ix<Nx; ix++) {
			exx = eps.xx({ix,iy,iz});
			eyy = eps.yy({ix,iy,iz});
			ezz = eps.zz({ix,iy,iz});
			exy = eps.xy({ix,iy,iz});
			exz = eps.xz({ix,iy,iz});
			eyz = eps.yz({ix,iy,iz});
			dst({ix,iy,0},dst_comp) =
				mu/(2*ezz*std::sqrt(eyy)) * (
					exy/2. * (
						d2_ovr_dY2(src, ix, iy, src_comp) -
						d2_ovr_dX2(src, ix, iy, src_comp) ) +
					(exx-ezz) * d2_ovr_dXdY(src, ix, iy, src_comp) ) -
				I*std::conj(mu)/(2*ezz) * (
					eyz*d_ovr_dX(src, ix, iy, src_comp) +
					exz*d_ovr_dY(src, ix, iy, src_comp) );
		}
	}
}
