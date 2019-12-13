#include "ParaxialPrimaryEvolutionOperator.h"

ParaxialPrimaryEvolutionOperator::ParaxialPrimaryEvolutionOperator(
		const PermittivityTensorField& eps,
		double wavelength, const RootSettings &settings) :
	BaseBPMOperator(eps, wavelength, settings.algorithm.bpm),
	d_ovr_dX_csr(Nx*Ny,Nx*Ny),
	d_ovr_dX_csc(Nx*Ny,Nx*Ny),
	d_ovr_dY_csr(Nx*Ny,Nx*Ny),
	d_ovr_dY_csc(Nx*Ny,Nx*Ny),
	D0_matrix(
		eps.mesh, {true, true, true, false},
		{false, true, true, true}, {true, false, false, true}),
	D1_matrix(
		eps.mesh, {false, true, true, true},
		{true, true, true, false}, {true, false, false, true}),
	D_matrix(eps.mesh),
	R_matrix(eps.mesh),
	N(Nx*Ny),
	N_woodbury_steps(settings.algorithm.bpm.n_woodbury_steps),
	tmp(eps.mesh),
	mu(std::complex<double>(0,delta_Z/2.)){

	std::vector<Eigen::Triplet<double> > triplet_list;
	int i;

	// First, we initialize the matrices d_ovr_dX_* 
	for(int iy=0; iy<Ny; iy++) {
		i = Nx*iy;
		triplet_list.push_back(Eigen::Triplet<double>(i, i+Nx-1, -0.5/delta_X));
		triplet_list.push_back(Eigen::Triplet<double>(i, i+1, 0.5/delta_X));
		for(int ix=1; ix<Nx-1; ix++) {
			i = ix+Nx*iy;
			triplet_list.push_back(Eigen::Triplet<double>(i, i-1, -0.5/delta_X));
			triplet_list.push_back(Eigen::Triplet<double>(i, i+1, 0.5/delta_X));
		}
		i = Nx-1+Nx*iy;
		triplet_list.push_back(Eigen::Triplet<double>(i, i-1, -0.5/delta_X));
		triplet_list.push_back(Eigen::Triplet<double>(i, Nx*iy, 0.5/delta_X));
	}
	d_ovr_dX_csr.setFromTriplets(triplet_list.begin(), triplet_list.end());
	d_ovr_dX_csc.setFromTriplets(triplet_list.begin(), triplet_list.end());

	// Second, we initialize the matrices d_ovr_dY_*.
	triplet_list.clear();
	for(int ix=0; ix<Nx; ix++) {
		i = ix;
		triplet_list.push_back(Eigen::Triplet<double>(i, i+Nx*(Ny-1), -0.5/delta_Y));
		triplet_list.push_back(Eigen::Triplet<double>(i, i+Nx, 0.5/delta_Y));
		for(int iy=1; iy<Ny-1; iy++) {
			i = ix+Nx*iy;
			triplet_list.push_back(Eigen::Triplet<double>(i, i-Nx, -0.5/delta_Y));
			triplet_list.push_back(Eigen::Triplet<double>(i, i+Nx, 0.5/delta_Y));
		}
		i = ix+Nx*(Ny-1);
		triplet_list.push_back(Eigen::Triplet<double>(i, i-Nx, -0.5/delta_Y));
		triplet_list.push_back(Eigen::Triplet<double>(i, ix, 0.5/delta_Y));
	}
	d_ovr_dY_csr.setFromTriplets(triplet_list.begin(), triplet_list.end());
	d_ovr_dY_csc.setFromTriplets(triplet_list.begin(), triplet_list.end());

	// Third, we initialize the matrix_D0.
	for(int iy=1; iy<Ny-1; iy++) {
		for(int ix=1; ix<Nx-1; ix++) {
			i = ix+Nx*iy;
			
			D0_matrix(DiffDY, Block00, {ix,iy}, myShift) = 1./std::pow(delta_Y,2.);
			D0_matrix(DiffDY, Block00, {ix,iy}, noShift) = -2./std::pow(delta_Y,2.);
			D0_matrix(DiffDY, Block00, {ix,iy}, pyShift) = 1./std::pow(delta_Y,2.);

			if(ix>=Nx_plate-1 && ix<Nx-Nx_plate+1) {
				D0_matrix(DiffDXDY, Block01, {ix,iy}, pxpyShift) = 0.25/(delta_X*delta_Y);
				D0_matrix(DiffDXDY, Block01, {ix,iy}, pxmyShift) = -0.25/(delta_X*delta_Y);
				D0_matrix(DiffDXDY, Block01, {ix,iy}, mxpyShift) = -0.25/(delta_X*delta_Y);
				D0_matrix(DiffDXDY, Block01, {ix,iy}, mxmyShift) = 0.25/(delta_X*delta_Y);

				D0_matrix(DiffDXDY, Block10, {ix,iy}, pxpyShift) = 0.25/(delta_X*delta_Y);
				D0_matrix(DiffDXDY, Block10, {ix,iy}, pxmyShift) = -0.25/(delta_X*delta_Y);
				D0_matrix(DiffDXDY, Block10, {ix,iy}, mxpyShift) = -0.25/(delta_X*delta_Y);
				D0_matrix(DiffDXDY, Block10, {ix,iy}, mxmyShift) = 0.25/(delta_X*delta_Y);
			}

			D0_matrix(DiffDX, Block11, {ix,iy}, mxShift) = 1./std::pow(delta_X,2.);
			D0_matrix(DiffDX, Block11, {ix,iy}, noShift) = -2./std::pow(delta_X,2.);
			D0_matrix(DiffDX, Block11, {ix,iy}, pxShift) = 1./std::pow(delta_X,2.);
		}
	}

	// Finally, we initialize the sample plate dofs of matrix_D1. LC dofs will be
	// initialized in the update_D1 method
	for(int iy=1; iy<Ny-1; iy++) {
		for(int ix=1; ix<Nx-1; ix++) {
			i = ix+Nx*iy;
			
			if(ix<Nx_plate-1 || ix>=Nx-Nx_plate+1) {
				D1_matrix(DiffDX, Block00, {ix,iy}, mxShift) = 1./std::pow(delta_X,2.);
				D1_matrix(DiffDX, Block00, {ix,iy}, noShift) = -2./std::pow(delta_X,2.);
				D1_matrix(DiffDX, Block00, {ix,iy}, pxShift) = 1./std::pow(delta_X,2.);

				D1_matrix(DiffDY, Block11, {ix,iy}, myShift) = 1./std::pow(delta_Y,2.);
				D1_matrix(DiffDY, Block11, {ix,iy}, noShift) = -2./std::pow(delta_Y,2.);
				D1_matrix(DiffDY, Block11, {ix,iy}, pyShift) = 1./std::pow(delta_Y,2.);
			}
		}
	}
}

int ParaxialPrimaryEvolutionOperator::apply(TransverseOpticalField &field) {

	// Forward y-step
	R_matrix.DY_op().block_vmult(Block00, tmp, 0, field, 0, false);
	R_matrix.DY_op().block_vmult(Block11, tmp, 1, field, 1, false);
	R_matrix.block_vmult(Block10, tmp, 1, field, 0, true);
	field.add_scaled_field(mu, tmp);


	// Backward x-step
	R_matrix.DX_op().shifted_block_vmult_inv(-mu, Block11, tmp, 1, field, 1);
	R_matrix.block_vmult(Block01, tmp, 0, tmp, 1);
	tmp.scale_and_add_block(mu, field, 0);
	R_matrix.DX_op().shifted_block_vmult_inv(-mu, Block00, field, 0, tmp, 0);
	field.copy_block(tmp, 1);

	tmp.copy_block(field, 0);
	for(int it=0; it<N_woodbury_steps; it++) {
		R_matrix.DXDY_op().block_vmult(Block00, tmp, 1, tmp, 0);
		tmp.scale_block(-mu, 1);
		R_matrix.DX_op().shifted_block_vmult_inv(-mu, Block00, tmp, 0, tmp, 1);
		field.add_block(tmp, 0);
	}


	// Forward x-step
	R_matrix.DX_op().block_vmult(Block00, tmp, 0, field, 0, false);
	R_matrix.DX_op().block_vmult(Block11, tmp, 1, field, 1, false);
	R_matrix.block_vmult(Block01, tmp, 0, field, 1, true);
	field.add_scaled_field(mu, tmp);


	// Backward y-step
	R_matrix.DY_op().shifted_block_vmult_inv(-mu, Block00, tmp, 0, field, 0);
	R_matrix.block_vmult(Block10, tmp, 1, tmp, 0);
	tmp.scale_and_add_block(mu, field, 1);
	R_matrix.DY_op().shifted_block_vmult_inv(-mu, Block11, field, 1, tmp, 1);
	field.copy_block(tmp, 0);

	tmp.copy_block(field, 1);
	for(int it=0; it<N_woodbury_steps; it++) {
		R_matrix.DXDY_op().block_vmult(Block11, tmp, 0, tmp, 1);
		tmp.scale_block(-mu, 0);
		R_matrix.DY_op().shifted_block_vmult_inv(-mu, Block11, tmp, 1, tmp, 0);
		field.add_block(tmp, 1);
	}

	return 2*N_woodbury_steps;
}

void ParaxialPrimaryEvolutionOperator::update() {

	// First we update the D1 and D  matrices with the values of permittivity for the current
	// transverse plane
	update_D1();
	update_D();

	// Then, we reinitialize R_matrix with I*W 
	typedef Eigen::SparseMatrix<double,Eigen::ColMajor>::InnerIterator ColIterator;
	typedef Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator RowIterator;

	for(int row=0; row<2*N; row++)
		for(auto diff_type : {DiffDX,DiffDY,DiffDXDY})
			for(auto dst_it : R_matrix.row_iterators(diff_type, row))
				dst_it.val = 0;

	#pragma omp parallel for
	for(int i=0; i<2*N; i++) {
		int mesh_idx = i%N;
		Index3D p({mesh_idx%Nx,mesh_idx/Nx,iz});
		Index3D p_pz({p.x,p.y,p.z+1});
		auto val = (i<N) ? 
			(eps.xz(p)/eps.zz(p)+eps.xz(p_pz)/eps.zz(p_pz))*std::complex<double>(0,0.25) :
			(eps.yz(p)/eps.zz(p)+eps.yz(p_pz)/eps.zz(p_pz))*std::complex<double>(0,0.25);

		// We add I*0.5*eps_az/eps_zz*\partial_b row-wise
		auto dst_row_it = R_matrix.row_iterators(DiffDX, i).begin();
		for(RowIterator src_it(d_ovr_dX_csr,mesh_idx); src_it; ++src_it) {
			while(dst_row_it->col<src_it.col())
				++dst_row_it;
			dst_row_it->val += src_it.value()*val;
		}
		dst_row_it = R_matrix.row_iterators(DiffDY, i).begin();
		for(RowIterator src_it(d_ovr_dY_csr,mesh_idx); src_it; ++src_it) {
			while(dst_row_it->col<src_it.col()+N)
				++dst_row_it;
			dst_row_it->val += src_it.value()*val;
		}

		// We add 0.5*\partial_a*eps_bz/eps_zz column-wise
		auto dst_col_it = R_matrix.column_iterators(DiffDX, i).begin();
		for(ColIterator src_it(d_ovr_dX_csc,mesh_idx); src_it; ++src_it) {
			while(dst_col_it->row<src_it.row())
				++dst_col_it;
			dst_col_it->val += src_it.value()*val;
		}
		dst_col_it = R_matrix.column_iterators(DiffDY, i).begin();
		for(ColIterator src_it(d_ovr_dY_csc,mesh_idx); src_it; ++src_it) {
			while(dst_col_it->row<src_it.row()+N)
				++dst_col_it;
			dst_col_it->val += src_it.value()*val;
		}
	}

	// Finally, we add (D*sqrt(K)^-1 + sqrt(K)^-1*D)/4 to R_matrix
	double exx_sqrt, exy_sqrt, eyy_sqrt, det, exx_sqrt_inv, exy_sqrt_inv, eyy_sqrt_inv;
	#pragma omp parallel for private(exx_sqrt,exy_sqrt,eyy_sqrt,det,exx_sqrt_inv, \
			exy_sqrt_inv,eyy_sqrt_inv)
	for(int mesh_idx=0; mesh_idx<N; mesh_idx++) {
		Index3D p({mesh_idx%Nx,mesh_idx/Nx,iz});
		Index3D p_pz({p.x,p.y,p.z+1});

		exx_sqrt = (eps.xx_sqrt(p)+eps.xx_sqrt(p_pz))/2;
		exy_sqrt = (eps.xy_sqrt(p)+eps.xy_sqrt(p_pz))/2;
		eyy_sqrt = (eps.yy_sqrt(p)+eps.yy_sqrt(p_pz))/2;
		det = exx_sqrt*eyy_sqrt - std::pow(exy_sqrt,2.);
		exx_sqrt_inv = eyy_sqrt/det;
		exy_sqrt_inv = -exy_sqrt/det;
		eyy_sqrt_inv = exx_sqrt/det;

		// We add 0.25*D*sqrt(K)^-1 to R_matrix column-wise
		for(auto diff_type : {DiffDX,DiffDY,DiffDXDY}) {
			// We add the left blocks (Block00 and Block10) of 0.25*D*sqrt(K)^-1 to R_matrix
			auto dst_it = R_matrix.column_iterators(diff_type, mesh_idx).begin();
			for(auto D_it : D_matrix.column_iterators(diff_type, mesh_idx)) {
				while(dst_it->row<D_it.row)
					++dst_it;
				dst_it->val += 0.25*D_it.val*exx_sqrt_inv;
			}
			dst_it = R_matrix.column_iterators(diff_type, mesh_idx).begin();
			for(auto D_it : D_matrix.column_iterators(diff_type, mesh_idx+N)) {
				while(dst_it->row<D_it.row)
					++dst_it;
				dst_it->val += 0.25*D_it.val*exy_sqrt_inv;
			}

			// We add the right blocks (Block01 and Block11) of D1*K to D_matrix
			dst_it = R_matrix.column_iterators(diff_type, mesh_idx+N).begin();
			for(auto D_it : D_matrix.column_iterators(diff_type, mesh_idx)) {
				while(dst_it->row<D_it.row)
					++dst_it;
				dst_it->val += 0.25*D_it.val*exy_sqrt_inv;
			}
			dst_it = R_matrix.column_iterators(diff_type, mesh_idx+N).begin();
			for(auto D_it : D_matrix.column_iterators(diff_type, mesh_idx+N)) {
				while(dst_it->row<D_it.row)
					++dst_it;
				dst_it->val += 0.25*D_it.val*eyy_sqrt_inv;
			}
		}

		// We add 0.25*sqrt(K)^-1*D to R_matrix row-wise
		for(auto diff_type : {DiffDX,DiffDY,DiffDXDY}) {
			// We add the upper blocks (Block00 and Block01) of 0.25*sqrt(K)^-1*D to R_matrix
			auto dst_it = R_matrix.row_iterators(diff_type, mesh_idx).begin();
			for(auto D_it : D_matrix.row_iterators(diff_type, mesh_idx)) {
				while(dst_it->col<D_it.col)
					++dst_it;
				dst_it->val += 0.25*exx_sqrt_inv*D_it.val;
			}
			dst_it = R_matrix.row_iterators(diff_type, mesh_idx).begin();
			for(auto D_it : D_matrix.row_iterators(diff_type, mesh_idx+N)) {
				while(dst_it->col<D_it.col)
					++dst_it;
				dst_it->val += 0.25*exy_sqrt_inv*D_it.val;
			}

			// We add the lower blocks (Block10 and Block11) of 0.25*sqrt(K)^-1*D to D_matrix
			dst_it = R_matrix.row_iterators(diff_type, mesh_idx+N).begin();
			for(auto D_it : D_matrix.row_iterators(diff_type, mesh_idx)) {
				while(dst_it->col<D_it.col)
					++dst_it;
				dst_it->val += 0.25*exy_sqrt_inv*D_it.val;
			}
			dst_it = R_matrix.row_iterators(diff_type, mesh_idx+N).begin();
			for(auto D_it : D_matrix.row_iterators(diff_type, mesh_idx+N)) {
				while(dst_it->col<D_it.col)
					++dst_it;
				dst_it->val += 0.25*eyy_sqrt_inv*D_it.val;
			}
		}
	}
}

void ParaxialPrimaryEvolutionOperator::update_D() {

	double exx_tr, exy_tr, eyy_tr;

	// We reinitialize D_matrix with D0_matrix
	#pragma omp parallel for
	for(int row=0; row<2*N; row++) {
		for(auto diff_type : {DiffDX,DiffDY,DiffDXDY}) {
			for(auto dst_it : D_matrix.row_iterators(diff_type, row))
				dst_it.val = 0;

			auto dst_it = D_matrix.row_iterators(diff_type, row).begin();
			for(auto D0_it : D0_matrix.row_iterators(diff_type, row)) {
				while(dst_it->col<D0_it.col)
					++dst_it;
				dst_it->val += D0_it.val;
			}
		}
	}

	// We add D1*K to D
	#pragma omp parallel for private(exx_tr, exy_tr, eyy_tr)
	for(int mesh_idx=0; mesh_idx<N; mesh_idx++) {
		Index3D p({mesh_idx%Nx,mesh_idx/Nx,iz});
		Index3D p_pz({p.x,p.y,p.z+1});

		exx_tr = (eps.xx_tr(p)+eps.xx_tr(p_pz))/2;
		exy_tr = (eps.xy_tr(p)+eps.xy_tr(p_pz))/2;
		eyy_tr = (eps.yy_tr(p)+eps.yy_tr(p_pz))/2;

		for(auto diff_type : {DiffDX,DiffDY,DiffDXDY}) {
			// We add the left blocks (Block00 and Block10) of D1*K to D_matrix
			auto dst_it = D_matrix.column_iterators(diff_type, mesh_idx).begin();
			for(auto D1_it : D1_matrix.column_iterators(diff_type, mesh_idx)) {
				while(dst_it->row<D1_it.row)
					++dst_it;
				dst_it->val += D1_it.val*exx_tr;
			}
			dst_it = D_matrix.column_iterators(diff_type, mesh_idx).begin();
			for(auto D1_it : D1_matrix.column_iterators(diff_type, mesh_idx+N)) {
				while(dst_it->row<D1_it.row)
					++dst_it;
				dst_it->val += D1_it.val*exy_tr;
			}

			// We add the right blocks (Block10 and Block11) of D1*K to D_matrix
			dst_it = D_matrix.column_iterators(diff_type, mesh_idx+N).begin();
			for(auto D1_it : D1_matrix.column_iterators(diff_type, mesh_idx)) {
				while(dst_it->row<D1_it.row)
					++dst_it;
				dst_it->val += D1_it.val*exy_tr;
			}
			dst_it = D_matrix.column_iterators(diff_type, mesh_idx+N).begin();
			for(auto D1_it : D1_matrix.column_iterators(diff_type, mesh_idx+N)) {
				while(dst_it->row<D1_it.row)
					++dst_it;
				dst_it->val += D1_it.val*eyy_tr;
			}
		}
	}
}

void ParaxialPrimaryEvolutionOperator::update_D1() {

	double ezz_inv, ezz_inv_px, ezz_inv_mx, ezz_inv_py, ezz_inv_my;
	#pragma omp parallel for private(ezz_inv,ezz_inv_px,ezz_inv_mx,ezz_inv_py,ezz_inv_my)
	for(int iy=1; iy<Ny-1; iy++) {
		for(int ix=Nx_plate-1; ix<Nx-Nx_plate+1; ix++) {
			Index3D p({ix,iy,iz});
			Index3D p_pz({ix,iy,iz+1});
			Index3D p_px({ix+1,iy,iz});
			Index3D p_px_pz({ix+1,iy,iz+1});
			Index3D p_mx({ix-1,iy,iz});
			Index3D p_mx_pz({ix-1,iy,iz+1});
			Index3D p_py({ix,iy+1,iz});
			Index3D p_py_pz({ix,iy+1,iz+1});
			Index3D p_my({ix,iy-1,iz});
			Index3D p_my_pz({ix,iy-1,iz+1});

			ezz_inv = (1./eps.zz(p)+1./eps.zz(p_pz))/2;
			ezz_inv_px = (1./eps.zz(p_px)+1./eps.zz(p_px_pz))/2;
			ezz_inv_mx = (1./eps.zz(p_mx)+1./eps.zz(p_mx_pz))/2;
			ezz_inv_py = (1./eps.zz(p_py)+1./eps.zz(p_py_pz))/2;
			ezz_inv_my = (1./eps.zz(p_my)+1./eps.zz(p_my_pz))/2;

			if(ix<=Nx_plate)
				D1_matrix(DiffDX, Block00, {ix,iy}, mxShift) = 
					ezz_inv_mx/std::pow(delta_X, 2.);
			else
				D1_matrix(DiffDX, Block00, {ix,iy}, mxShift) = 
					0.5*(ezz_inv_mx+ezz_inv)/std::pow(delta_X, 2.);
			if(ix>=Nx-Nx_plate-1)
				D1_matrix(DiffDX, Block00, {ix,iy}, pxShift) = 
					ezz_inv_px/std::pow(delta_X, 2.);
			else
				D1_matrix(DiffDX, Block00, {ix,iy}, pxShift) = 
					0.5*(ezz_inv_px+ezz_inv)/std::pow(delta_X, 2.);
			D1_matrix(DiffDX, Block00, {ix,iy}, noShift) = 
				- D1_matrix(DiffDX, Block00, {ix,iy}, mxShift)
				- D1_matrix(DiffDX, Block00, {ix,iy}, pxShift);

			D1_matrix(DiffDXDY, Block01, {ix,iy}, pxpyShift) =
				0.25*ezz_inv_px/(delta_X*delta_Y);
			D1_matrix(DiffDXDY, Block01, {ix,iy}, mxmyShift) =
				0.25*ezz_inv_mx/(delta_X*delta_Y);
			D1_matrix(DiffDXDY, Block01, {ix,iy}, mxpyShift) =
				- D1_matrix(DiffDXDY, Block01, {ix,iy}, mxmyShift);
			D1_matrix(DiffDXDY, Block01, {ix,iy}, pxmyShift) =
				- D1_matrix(DiffDXDY, Block01, {ix,iy}, pxpyShift);

			D1_matrix(DiffDXDY, Block10, {ix,iy}, pxpyShift) =
				0.25*ezz_inv_py/(delta_X*delta_Y);
			D1_matrix(DiffDXDY, Block10, {ix,iy}, mxmyShift) =
				0.25*ezz_inv_my/(delta_X*delta_Y);
			D1_matrix(DiffDXDY, Block10, {ix,iy}, mxpyShift) =
				- D1_matrix(DiffDXDY, Block10, {ix,iy}, pxpyShift);
			D1_matrix(DiffDXDY, Block10, {ix,iy}, pxmyShift) =
				- D1_matrix(DiffDXDY, Block10, {ix,iy}, mxmyShift);

			D1_matrix(DiffDY, Block11, {ix,iy}, pyShift) = 
				0.5*(ezz_inv_py+ezz_inv)/std::pow(delta_Y, 2.);
			D1_matrix(DiffDY, Block11, {ix,iy}, myShift) = 
				0.5*(ezz_inv_my+ezz_inv)/std::pow(delta_Y, 2.);
			D1_matrix(DiffDY, Block11, {ix,iy}, noShift) = 
				- D1_matrix(DiffDY, Block11, {ix,iy}, myShift)
				- D1_matrix(DiffDY, Block11, {ix,iy}, pyShift);
		}
	}
}
