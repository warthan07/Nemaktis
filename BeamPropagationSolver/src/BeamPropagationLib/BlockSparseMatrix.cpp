#include "BlockSparseMatrix.h"

template <typename T, int S>
BlockSparseMatrix<T,S>::BlockSparseMatrix(
		CartesianMesh mesh) :
	Nx(mesh.Nx),
	Ny(mesh.Ny),
	N(Nx*Ny),
	matrix(2*N,2*N) {

	col_iterator_arrays.resize(2*N);
	row_iterator_arrays.resize(2*N);
	for(int block_idx=0; block_idx<4; block_idx++)
		for(int shift_idx=0; shift_idx<9; shift_idx++)
			ptrs[block_idx][shift_idx].resize(N, NULL);
}

template <typename T, int S>
void BlockSparseMatrix<T,S>::update_pointers() {

	int xwrap = Nx-2;
	int ywrap = Nx*(Ny-2);

	int mesh_idx_row, flattened_shift, block_idx;
	Shift shift_idx;
	for(Eigen::Index k=0; k<matrix.outerSize(); k++) {
		for(typename Eigen::SparseMatrix<T,S>::InnerIterator it(matrix,k); it; ++it) {
			block_idx = (it.col()/N) + 2*(it.row()/N);
			flattened_shift = (it.col()%N) - (it.row()%N);

			if((it.row()%N)/Nx == (it.col()%N)/Nx) { // shifts along x
				if(flattened_shift==0)
					shift_idx = noShift;
				else if(flattened_shift==1 || flattened_shift==-xwrap)
					shift_idx = pxShift;
				else if(flattened_shift==-1 || flattened_shift==xwrap)
					shift_idx = mxShift;
			}
			else if((it.row()%N)%Nx == (it.col()%N)%Nx) { // shifts along y 
				if(flattened_shift==Nx || flattened_shift==-ywrap)
					shift_idx = pyShift;
				else if(flattened_shift==-Nx || flattened_shift==ywrap)
					shift_idx = myShift;
			}
			else { // shifts along both x and y
				if(flattened_shift==1+Nx || flattened_shift==1-ywrap ||
						flattened_shift==-xwrap+Nx || flattened_shift==-xwrap-ywrap)
					shift_idx = pxpyShift;
				else if(flattened_shift==-1+Nx || flattened_shift==-1-ywrap ||
						flattened_shift==xwrap+Nx || flattened_shift==xwrap-ywrap)
					shift_idx = mxpyShift;
				else if(flattened_shift==1-Nx || flattened_shift==1+ywrap ||
						flattened_shift==-xwrap-Nx || flattened_shift==-xwrap+ywrap)
					shift_idx = pxmyShift;
				else if(flattened_shift==-1-Nx || flattened_shift==-1+ywrap ||
						flattened_shift==xwrap-Nx || flattened_shift==xwrap+ywrap)
					shift_idx = mxmyShift;
			}
			
			ptrs[block_idx][shift_idx][it.row()%N] = &(it.valueRef());
			col_iterator_arrays[it.col()].push_back(
				ColIteratorEntry({it.row(),it.valueRef()}));
			row_iterator_arrays[it.row()].push_back(
				RowIteratorEntry({it.col(),it.valueRef()}));
		}
	}
}

template <typename T, int S>
BlockSparseMatrixDX<T,S>::BlockSparseMatrixDX(
		const CartesianMesh &mesh, const bool (&zero_blocks)[4]) :
	BlockSparseMatrix<T,S>(mesh) {

	std::vector<Eigen::Triplet<T> > triplet_list;
	int i, ipx, imx;

	auto Nx = this->Nx;
	auto Ny = this->Ny;
	auto N = this->N;

	// We only initialize the sparsity pattern if Nx>3 (for Nx<=3, we assume x-invariant
	// structures and fields)
	if(Nx>3) {
		for(int iy=0; iy<Ny; iy++) {
			for(int ix=0; ix<Nx-1; ix++) {
				i = ix+Nx*iy;
				ipx = (ix==Nx-2) ? i+2-Nx : i+1;
				imx = (ix==0) ? i+Nx-2 : i-1;
	
				if(!zero_blocks[Block00]) {
					triplet_list.push_back(Eigen::Triplet<T>(i, imx, 0));
					triplet_list.push_back(Eigen::Triplet<T>(i, i, 0));
					triplet_list.push_back(Eigen::Triplet<T>(i, ipx, 0));
				}
	
				if(!zero_blocks[Block01]) {
					triplet_list.push_back(Eigen::Triplet<T>(i, N+imx, 0));
					triplet_list.push_back(Eigen::Triplet<T>(i, N+i, 0));
					triplet_list.push_back(Eigen::Triplet<T>(i, N+ipx, 0));
				}
	
				if(!zero_blocks[Block10]) {
					triplet_list.push_back(Eigen::Triplet<T>(N+i, imx, 0));
					triplet_list.push_back(Eigen::Triplet<T>(N+i, i, 0));
					triplet_list.push_back(Eigen::Triplet<T>(N+i, ipx, 0));
				}
	
				if(!zero_blocks[Block11]) {
					triplet_list.push_back(Eigen::Triplet<T>(N+i, N+imx, 0));
					triplet_list.push_back(Eigen::Triplet<T>(N+i, N+i, 0));
					triplet_list.push_back(Eigen::Triplet<T>(N+i, N+ipx, 0));
				}
			}
		}
	}
	this->matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
	this->update_pointers();
}

template <typename T, int S>
void BlockSparseMatrixDX<T,S>::block_vmult(
		const Block &block_idx,
		TransverseOpticalField &dst, int dst_comp,
		TransverseOpticalField &src, int src_comp,
		bool add) const {

	if(add) {
		if(this->Nx>3) {
			dst().segment(dst_comp*this->N,this->N) +=
				this->matrix.block(
					(block_idx/2)*this->N,(block_idx%2)*this->N,this->N,this->N)*
						src().segment(src_comp*this->N,this->N);
			#pragma omp parallel for
			for(int iy=0; iy<this->Ny; iy++) {
				// int ix = 0;
				// dst({ix,iy},dst_comp) +=
					// (*this)(block_idx,{ix,iy},mxShift) * src({ix+this->Nx-2,iy},src_comp) +
					// (*this)(block_idx,{ix,iy},noShift) * src({ix,iy},src_comp) +
					// (*this)(block_idx,{ix,iy},pxShift) * src({ix+1,iy},src_comp);
				// for(ix=1; ix<this->Nx-1; ix++) {
					// dst({ix,iy},dst_comp) +=
						// (*this)(block_idx,{ix,iy},mxShift) * src({ix-1,iy},src_comp) +
						// (*this)(block_idx,{ix,iy},noShift) * src({ix,iy},src_comp) +
						// (*this)(block_idx,{ix,iy},pxShift) * src({ix+1,iy},src_comp);
				// }
				// dst({ix,iy},dst_comp) = dst({0,iy},dst_comp);
				dst({this->Nx-1,iy},dst_comp) = dst({0,iy},dst_comp);
			}
		}
	}
	else {
		if(this->Nx>3) {
			dst().segment(dst_comp*this->N,this->N) =
				this->matrix.block(
					(block_idx/2)*this->N,(block_idx%2)*this->N,this->N,this->N)*
						src().segment(src_comp*this->N,this->N);
			#pragma omp parallel for
			for(int iy=0; iy<this->Ny; iy++) {
				// int ix = 0;
				// dst({ix,iy},dst_comp) =
					// (*this)(block_idx,{ix,iy},mxShift) * src({ix+this->Nx-2,iy},src_comp) +
					// (*this)(block_idx,{ix,iy},noShift) * src({ix,iy},src_comp) +
					// (*this)(block_idx,{ix,iy},pxShift) * src({ix+1,iy},src_comp);
				// for(ix=1; ix<this->Nx-1; ix++) {
					// dst({ix,iy},dst_comp) =
						// (*this)(block_idx,{ix,iy},mxShift) * src({ix-1,iy},src_comp) +
						// (*this)(block_idx,{ix,iy},noShift) * src({ix,iy},src_comp) +
						// (*this)(block_idx,{ix,iy},pxShift) * src({ix+1,iy},src_comp);
				// }
				// dst({ix,iy},dst_comp) = dst({0,iy},dst_comp);
				dst({this->Nx-1,iy},dst_comp) = dst({0,iy},dst_comp);
			}
		}
		else
			dst.set_block_to_zero(dst_comp);
	}
}

template <typename T, int S>
void BlockSparseMatrixDX<T,S>::shifted_block_vmult(
		std::complex<double> scaling_factor, const Block &block_idx,
		TransverseOpticalField &dst, int dst_comp,
		TransverseOpticalField &src, int src_comp,
		bool add) const {

	std::complex<double> val;

	if(add) {
		if(this->Nx>3) {
			dst().segment(dst_comp*this->N,this->N) +=
				src().segment(src_comp*this->N,this->N) +
				scaling_factor*this->matrix.block(
					(block_idx/2)*this->N,(block_idx%2)*this->N,this->N,this->N)*
						src().segment(src_comp*this->N,this->N);
			#pragma omp parallel for
			for(int iy=0; iy<this->Ny; iy++) {
				// int ix = 0;
				// val =
					// (*this)(block_idx,{ix,iy},mxShift) * src({ix+this->Nx-2,iy},src_comp) +
					// (*this)(block_idx,{ix,iy},noShift) * src({ix,iy},src_comp) +
					// (*this)(block_idx,{ix,iy},pxShift) * src({ix+1,iy},src_comp);
				// dst({ix,iy},dst_comp) +=
					// src({ix,iy},src_comp) + scaling_factor*val;
				// for(ix=1; ix<this->Nx-1; ix++) {
					// val =
						// (*this)(block_idx,{ix,iy},mxShift) * src({ix-1,iy},src_comp) +
						// (*this)(block_idx,{ix,iy},noShift) * src({ix,iy},src_comp) +
						// (*this)(block_idx,{ix,iy},pxShift) * src({ix+1,iy},src_comp);
					// dst({ix,iy},dst_comp) +=
						// src({ix,iy},src_comp) + scaling_factor*val;
				// }
				// dst({ix,iy},dst_comp) = dst({0,iy},dst_comp);
				dst({this->Nx-1,iy},dst_comp) = dst({0,iy},dst_comp);
			}
		}
		else
			dst.add_block(dst_comp, src, src_comp);
	}
	else {
		if(this->Nx>3) {
			dst().segment(dst_comp*this->N,this->N) =
				src().segment(src_comp*this->N,this->N) +
				scaling_factor*this->matrix.block(
					(block_idx/2)*this->N,(block_idx%2)*this->N,this->N,this->N)*
						src().segment(src_comp*this->N,this->N);
			#pragma omp parallel for
			for(int iy=0; iy<this->Ny; iy++) {
				// int ix = 0;
				// val =
					// (*this)(block_idx,{ix,iy},mxShift) * src({ix+this->Nx-2,iy},src_comp) +
					// (*this)(block_idx,{ix,iy},noShift) * src({ix,iy},src_comp) +
					// (*this)(block_idx,{ix,iy},pxShift) * src({ix+1,iy},src_comp);
				// dst({ix,iy},dst_comp) =
					// src({ix,iy},src_comp) + scaling_factor*val;
				// for(ix=1; ix<this->Nx-1; ix++) {
					// val =
						// (*this)(block_idx,{ix,iy},mxShift) * src({ix-1,iy},src_comp) +
						// (*this)(block_idx,{ix,iy},noShift) * src({ix,iy},src_comp) +
						// (*this)(block_idx,{ix,iy},pxShift) * src({ix+1,iy},src_comp);
					// dst({ix,iy},dst_comp) =
						// src({ix,iy},src_comp) + scaling_factor*val;
				// }
				// dst({ix,iy},dst_comp) = dst({0,iy},dst_comp);
				dst({this->Nx-1,iy},dst_comp) = dst({0,iy},dst_comp);
			}
		}
		else
			dst.copy_block(dst_comp, src, src_comp);
	}
}

template <typename T, int S>
void BlockSparseMatrixDX<T,S>::shifted_block_vmult_inv(
		std::complex<double> scaling_factor, const Block &block_idx,
		TransverseOpticalField &dst, int dst_comp,
		const TransverseOpticalField &src, int src_comp) const {

	auto Nx = this->Nx;
	auto Ny = this->Ny;

	if(Nx<=3) {
		dst.copy_block(dst_comp, src, src_comp);
		return;
	}

	#pragma omp parallel for
	for(int iy=0; iy<Ny; iy++) {
		std::vector<std::complex<double> > gamma(Nx-2);
		std::vector<std::complex<double> > corr(Nx-1);

		std::complex<double> beta = 2.*diag_mid(scaling_factor, block_idx, {0,iy});
		dst({0,iy},dst_comp) = src({0,iy},src_comp) / beta;
		corr[0] = -0.5;
		
		for(int ix=1; ix<Nx-2; ix++) {
			gamma[ix-1] = diag_sup(scaling_factor, block_idx, {ix-1,iy}) / beta;
			beta =
				diag_mid(scaling_factor,block_idx,{ix,iy}) -
				diag_inf(scaling_factor,block_idx,{ix,iy})*gamma[ix-1];

			dst({ix,iy},dst_comp) =
				( src({ix,iy},src_comp) 
				- diag_inf(scaling_factor,block_idx,{ix,iy})*dst({ix-1,iy},dst_comp) ) / beta;
			corr[ix] = - diag_inf(scaling_factor,block_idx,{ix,iy})*corr[ix-1] / beta;
		}

		gamma[Nx-3] = diag_sup(scaling_factor, block_idx, {Nx-3,iy}) / beta;
		beta =
			diag_mid(scaling_factor, block_idx, {Nx-2,iy}) +
			diag_sup(scaling_factor, block_idx, {Nx-2,iy})
				* diag_inf(scaling_factor, block_idx, {0,iy}) 
				/ diag_mid(scaling_factor, block_idx, {0,iy})
			- diag_inf(scaling_factor, block_idx, {Nx-2,iy})*gamma[Nx-3];

		dst({Nx-2,iy},dst_comp) =
			( src({Nx-2,iy},src_comp)
			- diag_inf(scaling_factor, block_idx, {Nx-2,iy})*dst({Nx-3,iy},dst_comp) ) / beta;
		corr[Nx-2] = 
			( diag_sup(scaling_factor, block_idx, {Nx-2,iy})
			- diag_inf(scaling_factor, block_idx, {Nx-2,iy})*corr[Nx-3] ) / beta;

		for(int ix=Nx-2; ix>0; ix--) {
			dst({ix-1,iy},dst_comp) -= gamma[ix-1]*dst({ix,iy},dst_comp);
			corr[ix-1] -= gamma[ix-1]*corr[ix];
		}

		std::complex<double> frac =
			( dst({0,iy},dst_comp) - diag_inf(scaling_factor, block_idx, {0,iy})
				/ diag_mid(scaling_factor, block_idx, {0,iy})*dst({Nx-2,iy},dst_comp) ) /
			( 1. + corr[0] - diag_inf(scaling_factor, block_idx, {0,iy})
			  	/ diag_mid(scaling_factor, block_idx, {0,iy})*corr[Nx-2] );
		for(int ix=0; ix<Nx-1; ix++)
			dst({ix,iy},dst_comp) -= frac*corr[ix];
		dst({Nx-1,iy},dst_comp) = dst({0,iy},dst_comp);
	}
}

template <typename T, int S>
BlockSparseMatrixDY<T,S>::BlockSparseMatrixDY(
		const CartesianMesh &mesh, const bool (&zero_blocks)[4]) :
	BlockSparseMatrix<T,S>(mesh) {

	std::vector<Eigen::Triplet<T> > triplet_list;
	int i, imy, ipy;

	auto Nx = this->Nx;
	auto Ny = this->Ny;
	auto N = this->N;

	// We only initialize the sparsity pattern if Ny>3 (for Ny<=3, we assume y-invariant
	// structures and fields)
	if(Ny>3) {
		for(int ix=0; ix<Nx; ix++) {
			for(int iy=0; iy<Ny-1; iy++) {
				i = ix+Nx*iy;
				ipy = (iy==Ny-2) ? i+Nx*(2-Ny) : i+Nx;
				imy = (iy==0) ? i+Nx*(Ny-2) : i-Nx;
	
				if(!zero_blocks[Block00]) {
					triplet_list.push_back(Eigen::Triplet<T>(i, imy, 0));
					triplet_list.push_back(Eigen::Triplet<T>(i, i, 0));
					triplet_list.push_back(Eigen::Triplet<T>(i, ipy, 0));
				}
	
				if(!zero_blocks[Block01]) {
					triplet_list.push_back(Eigen::Triplet<T>(i, N+imy, 0));
					triplet_list.push_back(Eigen::Triplet<T>(i, N+i, 0));
					triplet_list.push_back(Eigen::Triplet<T>(i, N+ipy, 0));
				}
	
				if(!zero_blocks[Block10]) {
					triplet_list.push_back(Eigen::Triplet<T>(N+i, imy, 0));
					triplet_list.push_back(Eigen::Triplet<T>(N+i, i, 0));
					triplet_list.push_back(Eigen::Triplet<T>(N+i, ipy, 0));
				}
	
				if(!zero_blocks[Block11]) {
					triplet_list.push_back(Eigen::Triplet<T>(N+i, N+imy, 0));
					triplet_list.push_back(Eigen::Triplet<T>(N+i, N+i, 0));
					triplet_list.push_back(Eigen::Triplet<T>(N+i, N+ipy, 0));
				}
			}
		}
	}
	this->matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
	this->update_pointers();
}

template <typename T, int S>
void BlockSparseMatrixDY<T,S>::block_vmult(
		const Block &block_idx,
		TransverseOpticalField &dst, int dst_comp,
		TransverseOpticalField &src, int src_comp,
		bool add) const {

	if(add) {
		if(this->Ny>3) {
			dst().segment(dst_comp*this->N,this->N) +=
				this->matrix.block(
					(block_idx/2)*this->N,(block_idx%2)*this->N,this->N,this->N)*
						src().segment(src_comp*this->N,this->N);
			#pragma omp parallel for
			for(int ix=0; ix<this->Nx; ix++) {
				// int iy = 0;
				// dst({ix,iy},dst_comp) +=
					// (*this)(block_idx,{ix,iy},myShift) * src({ix,iy+this->Ny-2},src_comp) +
					// (*this)(block_idx,{ix,iy},noShift) * src({ix,iy},src_comp) +
					// (*this)(block_idx,{ix,iy},pyShift) * src({ix,iy+1},src_comp);
				// for(iy=1; iy<this->Ny-1; iy++) {
					// dst({ix,iy},dst_comp) +=
						// (*this)(block_idx,{ix,iy},myShift) * src({ix,iy-1},src_comp) +
						// (*this)(block_idx,{ix,iy},noShift) * src({ix,iy},src_comp) +
						// (*this)(block_idx,{ix,iy},pyShift) * src({ix,iy+1},src_comp);
				// }
				// dst({ix,iy},dst_comp) = dst({ix,0},dst_comp);
				dst({ix,this->Ny-1},dst_comp) = dst({ix,0},dst_comp);
			}
		}
	}
	else {
		if(this->Ny>3) {
			dst().segment(dst_comp*this->N,this->N) =
				this->matrix.block(
					(block_idx/2)*this->N,(block_idx%2)*this->N,this->N,this->N)*
						src().segment(src_comp*this->N,this->N);
			#pragma omp parallel for
			for(int ix=0; ix<this->Nx; ix++) {
				// int iy = 0;
				// dst({ix,iy},dst_comp) =
					// (*this)(block_idx,{ix,iy},myShift) * src({ix,iy+this->Ny-2},src_comp) +
					// (*this)(block_idx,{ix,iy},noShift) * src({ix,iy},src_comp) +
					// (*this)(block_idx,{ix,iy},pyShift) * src({ix,iy+1},src_comp);
				// for(iy=1; iy<this->Ny-1; iy++) {
					// dst({ix,iy},dst_comp) =
						// (*this)(block_idx,{ix,iy},myShift) * src({ix,iy-1},src_comp) +
						// (*this)(block_idx,{ix,iy},noShift) * src({ix,iy},src_comp) +
						// (*this)(block_idx,{ix,iy},pyShift) * src({ix,iy+1},src_comp);
				// }
				// dst({ix,iy},dst_comp) = dst({ix,0},dst_comp);
				dst({ix,this->Ny-1},dst_comp) = dst({ix,0},dst_comp);
			}
		}
		else
			dst.set_block_to_zero(dst_comp);
	}
}

template <typename T, int S>
void BlockSparseMatrixDY<T,S>::shifted_block_vmult(
		std::complex<double> scaling_factor, const Block &block_idx,
		TransverseOpticalField &dst, int dst_comp,
		TransverseOpticalField &src, int src_comp,
		bool add) const {

	std::complex<double> val;

	if(add) {
		if(this->Ny>3) {
			dst().segment(dst_comp*this->N,this->N) +=
				src().segment(src_comp*this->N,this->N) +
				scaling_factor*this->matrix.block(
					(block_idx/2)*this->N,(block_idx%2)*this->N,this->N,this->N)*
						src().segment(src_comp*this->N,this->N);
			#pragma omp parallel for
			for(int ix=0; ix<this->Nx; ix++) {
				// int iy = 0;
				// val =
					// (*this)(block_idx,{ix,iy},myShift) * src({ix,iy+this->Ny-2},src_comp) +
					// (*this)(block_idx,{ix,iy},noShift) * src({ix,iy},src_comp) +
					// (*this)(block_idx,{ix,iy},pyShift) * src({ix,iy+1},src_comp);
				// dst({ix,iy},dst_comp) +=
					// src({ix,iy},src_comp) + scaling_factor*val;
				// for(int iy=1; iy<this->Ny-1; iy++) {
					// val =
						// (*this)(block_idx,{ix,iy},myShift) * src({ix,iy-1},src_comp) +
						// (*this)(block_idx,{ix,iy},noShift) * src({ix,iy},src_comp) +
						// (*this)(block_idx,{ix,iy},pyShift) * src({ix,iy+1},src_comp);
					// dst({ix,iy},dst_comp) +=
						// src({ix,iy},src_comp) + scaling_factor*val;
				// }
				// dst({ix,iy},dst_comp) = dst({ix,0},dst_comp);
				dst({ix,this->Ny-1},dst_comp) = dst({ix,0},dst_comp);
			}
		}
		else
			dst.add_block(dst_comp, src, src_comp);
	}
	else {
		if(this->Ny>3) {
			dst().segment(dst_comp*this->N,this->N) =
				src().segment(src_comp*this->N,this->N) +
				scaling_factor*this->matrix.block(
					(block_idx/2)*this->N,(block_idx%2)*this->N,this->N,this->N)*
						src().segment(src_comp*this->N,this->N);
			#pragma omp parallel for
			for(int ix=0; ix<this->Nx; ix++) {
				// int iy = 0;
				// val =
					// (*this)(block_idx,{ix,iy},myShift) * src({ix,iy+this->Ny-2},src_comp) +
					// (*this)(block_idx,{ix,iy},noShift) * src({ix,iy},src_comp) +
					// (*this)(block_idx,{ix,iy},pyShift) * src({ix,iy+1},src_comp);
				// dst({ix,iy},dst_comp) =
					// src({ix,iy},src_comp) + scaling_factor*val;
				// for(iy=1; iy<this->Ny-1; iy++) {
					// val =
						// (*this)(block_idx,{ix,iy},myShift) * src({ix,iy-1},src_comp) +
						// (*this)(block_idx,{ix,iy},noShift) * src({ix,iy},src_comp) +
						// (*this)(block_idx,{ix,iy},pyShift) * src({ix,iy+1},src_comp);
					// dst({ix,iy},dst_comp) =
						// src({ix,iy},src_comp) + scaling_factor*val;
				// }
				// dst({ix,iy},dst_comp) = dst({ix,0},dst_comp);
				dst({ix,this->Ny-1},dst_comp) = dst({ix,0},dst_comp);
			}
		}
		else 
			dst.copy_block(dst_comp, src, src_comp);
	}
}

template <typename T, int S>
void BlockSparseMatrixDY<T,S>::shifted_block_vmult_inv(
		std::complex<double> scaling_factor, const Block &block_idx,
		TransverseOpticalField &dst, int dst_comp,
		const TransverseOpticalField &src, int src_comp) const {

	auto Nx = this->Nx;
	auto Ny = this->Ny;

	if(Ny<=3) {
		dst.copy_block(dst_comp, src, src_comp);
		return;
	}

	#pragma omp parallel for
	for(int ix=0; ix<Nx; ix++) {
		std::vector<std::complex<double> > gamma(Ny-2);
		std::vector<std::complex<double> > corr(Ny-1);

		std::complex<double> beta = 2.*diag_mid(scaling_factor, block_idx, {ix,0});
		dst({ix,0},dst_comp) = src({ix,0},src_comp) / beta;
		corr[0] = -0.5;

		for(int iy=1; iy<Ny-2; iy++) {
			gamma[iy-1] = diag_sup(scaling_factor, block_idx, {ix,iy-1}) / beta;
			beta =
				diag_mid(scaling_factor, block_idx, {ix,iy}) -
				diag_inf(scaling_factor, block_idx, {ix,iy})*gamma[iy-1];

			dst({ix,iy},dst_comp) =
				( src({ix,iy},src_comp) -
				diag_inf(scaling_factor, block_idx, {ix,iy})*dst({ix,iy-1},dst_comp)) / beta;
			corr[iy] = - diag_inf(scaling_factor, block_idx, {ix,iy})*corr[iy-1] / beta;
		}

		gamma[Ny-3] = diag_sup(scaling_factor, block_idx, {ix,Ny-3}) / beta;
		beta =
			diag_mid(scaling_factor, block_idx, {ix,Ny-2}) +
			diag_sup(scaling_factor, block_idx, {ix,Ny-2})
				* diag_inf(scaling_factor, block_idx, {ix,0})
				/ diag_mid(scaling_factor, block_idx, {ix,0})
			- diag_inf(scaling_factor, block_idx, {ix,Ny-2})*gamma[Ny-3];

		dst({ix,Ny-2},dst_comp) =
			( src({ix,Ny-2},src_comp) 
			- diag_inf(scaling_factor, block_idx, {ix,Ny-2})*dst({ix,Ny-3},dst_comp) ) / beta;
		corr[Ny-2] =
			( diag_sup(scaling_factor, block_idx, {ix,Ny-2})
			- diag_inf(scaling_factor, block_idx, {ix,Ny-2})*corr[Ny-3] ) / beta;

		for(int iy=Ny-2; iy>0; iy--) {
			dst({ix,iy-1},dst_comp) -= gamma[iy-1]*dst({ix,iy},dst_comp);
			corr[iy-1] -= gamma[iy-1]*corr[iy];
		}

		std::complex<double> frac =
			( dst({ix,0},dst_comp)
			- diag_inf(scaling_factor, block_idx, {ix,0})
				/ diag_mid(scaling_factor, block_idx, {ix,0})*dst({ix,Ny-2},dst_comp) ) /
			( 1. + corr[0] - diag_inf(scaling_factor, block_idx, {ix,0})
				/ diag_mid(scaling_factor, block_idx, {ix,0})*corr[Ny-2] );
		for(int iy=0; iy<Ny-1; iy++)
			dst({ix,iy},dst_comp) -= frac*corr[iy];
		dst({ix,Ny-1},dst_comp) = dst({ix,0},dst_comp);
	}
}

template <typename T, int S>
BlockSparseMatrixDXDY<T,S>::BlockSparseMatrixDXDY(
		const CartesianMesh &mesh, const bool (&zero_blocks)[4]) :
	BlockSparseMatrix<T,S>(mesh) {

	std::vector<Eigen::Triplet<T> > triplet_list;
	int px, py, mx, my, i;

	auto Nx = this->Nx;
	auto Ny = this->Ny;
	auto N = this->N;

	// We only initialize the sparsity pattern if Nx>3 and Ny>3 (else, we assume either
	// x-invariant or y-invariant structures and fields)
	if(Nx>3 && Ny>3) {
		for(int iy=0; iy<this->Ny-1; iy++) {
			for(int ix=0; ix<this->Nx-1; ix++) {
				i = ix+this->Nx*iy;
				px = (ix==Nx-2) ? 2-Nx : 1;
				mx = (ix==0) ? Nx-2 : -1;
				py = (iy==Ny-2) ? Nx*(2-Ny) : Nx;
				my = (iy==0) ? Nx*(Ny-2) : -Nx;
	
				if(!zero_blocks[Block00]) {
					triplet_list.push_back(Eigen::Triplet<T>(i, i+px+py, 0));
					triplet_list.push_back(Eigen::Triplet<T>(i, i+mx+py, 0));
					triplet_list.push_back(Eigen::Triplet<T>(i, i+px+my, 0));
					triplet_list.push_back(Eigen::Triplet<T>(i, i+mx+my, 0));
				}
	
				if(!zero_blocks[Block01]) {
					triplet_list.push_back(Eigen::Triplet<T>(i, N+i+px+py, 0));
					triplet_list.push_back(Eigen::Triplet<T>(i, N+i+mx+py, 0));
					triplet_list.push_back(Eigen::Triplet<T>(i, N+i+px+my, 0));
					triplet_list.push_back(Eigen::Triplet<T>(i, N+i+mx+my, 0));
				}
	
				if(!zero_blocks[Block10]) {
					triplet_list.push_back(Eigen::Triplet<T>(N+i, i+px+py, 0));
					triplet_list.push_back(Eigen::Triplet<T>(N+i, i+mx+py, 0));
					triplet_list.push_back(Eigen::Triplet<T>(N+i, i+px+my, 0));
					triplet_list.push_back(Eigen::Triplet<T>(N+i, i+mx+my, 0));
				}
	
				if(!zero_blocks[Block11]) {
					triplet_list.push_back(Eigen::Triplet<T>(N+i, N+i+px+py, 0));
					triplet_list.push_back(Eigen::Triplet<T>(N+i, N+i+mx+py, 0));
					triplet_list.push_back(Eigen::Triplet<T>(N+i, N+i+px+my, 0));
					triplet_list.push_back(Eigen::Triplet<T>(N+i, N+i+mx+my, 0));
				}
			}
		}
	}
	this->matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
	this->update_pointers();
}

template <typename T, int S>
void BlockSparseMatrixDXDY<T,S>::block_vmult(
		const Block &block_idx,
		TransverseOpticalField &dst, int dst_comp,
		TransverseOpticalField &src, int src_comp,
		bool add) const {

	if(add) {
		if(this->Nx>3 && this->Ny>3) {
			dst().segment(dst_comp*this->N,this->N) +=
				this->matrix.block(
					(block_idx/2)*this->N,(block_idx%2)*this->N,this->N,this->N)*
						src().segment(src_comp*this->N,this->N);
			#pragma omp parallel for
			for(int iy=1; iy<this->Ny-1; iy++) {
				// int ix = 0;
				// dst({ix,iy},dst_comp) +=
					// (*this)(block_idx,{ix,iy},pxpyShift) * src({ix+1,iy+1},src_comp) +
					// (*this)(block_idx,{ix,iy},pxmyShift) * src({ix+1,iy-1},src_comp) +
					// (*this)(block_idx,{ix,iy},mxpyShift) * src({this->Nx-2,iy+1},src_comp) +
					// (*this)(block_idx,{ix,iy},mxmyShift) * src({this->Nx-2,iy-1},src_comp);
				// for(ix=1; ix<this->Nx-1; ix++) {
					// dst({ix,iy},dst_comp) +=
						// (*this)(block_idx,{ix,iy},pxpyShift) * src({ix+1,iy+1},src_comp) +
						// (*this)(block_idx,{ix,iy},pxmyShift) * src({ix+1,iy-1},src_comp) +
						// (*this)(block_idx,{ix,iy},mxpyShift) * src({ix-1,iy+1},src_comp) +
						// (*this)(block_idx,{ix,iy},mxmyShift) * src({ix-1,iy-1},src_comp);
				// }
				// dst({ix,iy},dst_comp) = dst({0,iy}, dst_comp);
				dst({this->Nx-1,iy},dst_comp) = dst({0,iy}, dst_comp);
			}
			// int ix = 0;
			// dst({ix,0},dst_comp) +=
				// (*this)(block_idx,{ix,0},pxpyShift) * src({ix+1,1},src_comp) +
				// (*this)(block_idx,{ix,0},pxmyShift) * src({ix+1,this->Ny-2},src_comp) +
				// (*this)(block_idx,{ix,0},mxpyShift) * src({this->Nx-2,1},src_comp) +
				// (*this)(block_idx,{ix,0},mxmyShift) * src({this->Nx-2,this->Ny-2},src_comp);
			// for(ix=1; ix<this->Nx-1; ix++) {
				// dst({ix,0},dst_comp) +=
					// (*this)(block_idx,{ix,0},pxpyShift) * src({ix+1,1},src_comp) +
					// (*this)(block_idx,{ix,0},pxmyShift) * src({ix+1,+this->Ny-2},src_comp) +
					// (*this)(block_idx,{ix,0},mxpyShift) * src({ix-1,1},src_comp) +
					// (*this)(block_idx,{ix,0},mxmyShift) * src({ix-1,+this->Ny-2},src_comp);
			// }
			// dst({ix,0},dst_comp) = dst({0,0}, dst_comp);
			// for(ix=0; ix<this->Nx; ix++)
				// dst({ix,this->Ny-1},dst_comp) = dst({ix,0}, dst_comp);
			dst({this->Nx-1,0},dst_comp) = dst({0,0}, dst_comp);
			#pragma omp parallel for
			for(int ix=0; ix<this->Nx; ix++)
				dst({ix,this->Ny-1},dst_comp) = dst({ix,0}, dst_comp);
		}
	}
	else {
		if(this->Nx>3 && this->Ny>3) {
			dst().segment(dst_comp*this->N,this->N) =
				this->matrix.block(
					(block_idx/2)*this->N,(block_idx%2)*this->N,this->N,this->N)*
						src().segment(src_comp*this->N,this->N);
			#pragma omp parallel for
			for(int iy=1; iy<this->Ny-1; iy++) {
				// int ix = 0;
				// dst({ix,iy},dst_comp) =
					// (*this)(block_idx,{ix,iy},pxpyShift) * src({ix+1,iy+1},src_comp) +
					// (*this)(block_idx,{ix,iy},pxmyShift) * src({ix+1,iy-1},src_comp) +
					// (*this)(block_idx,{ix,iy},mxpyShift) * src({this->Nx-2,iy+1},src_comp) +
					// (*this)(block_idx,{ix,iy},mxmyShift) * src({this->Nx-2,iy-1},src_comp);
				// for(ix=1; ix<this->Nx-1; ix++) {
					// dst({ix,iy},dst_comp) =
						// (*this)(block_idx,{ix,iy},pxpyShift) * src({ix+1,iy+1},src_comp) +
						// (*this)(block_idx,{ix,iy},pxmyShift) * src({ix+1,iy-1},src_comp) +
						// (*this)(block_idx,{ix,iy},mxpyShift) * src({ix-1,iy+1},src_comp) +
						// (*this)(block_idx,{ix,iy},mxmyShift) * src({ix-1,iy-1},src_comp);
				// }
				// dst({ix,iy},dst_comp) = dst({0,iy}, dst_comp);
				dst({this->Nx-1,iy},dst_comp) = dst({0,iy}, dst_comp);
			}
			// int ix = 0;
			// dst({ix,0},dst_comp) =
				// (*this)(block_idx,{ix,0},pxpyShift) * src({ix+1,1},src_comp) +
				// (*this)(block_idx,{ix,0},pxmyShift) * src({ix+1,this->Ny-2},src_comp) +
				// (*this)(block_idx,{ix,0},mxpyShift) * src({this->Nx-2,1},src_comp) +
				// (*this)(block_idx,{ix,0},mxmyShift) * src({this->Nx-2,this->Ny-2},src_comp);
			// for(ix=1; ix<this->Nx-1; ix++) {
				// dst({ix,0},dst_comp) =
					// (*this)(block_idx,{ix,0},pxpyShift) * src({ix+1,1},src_comp) +
					// (*this)(block_idx,{ix,0},pxmyShift) * src({ix+1,+this->Ny-2},src_comp) +
					// (*this)(block_idx,{ix,0},mxpyShift) * src({ix-1,1},src_comp) +
					// (*this)(block_idx,{ix,0},mxmyShift) * src({ix-1,+this->Ny-2},src_comp);
			// }
			// dst({ix,0},dst_comp) = dst({0,0}, dst_comp);
			// for(ix=0; ix<this->Nx; ix++)
				// dst({ix,this->Ny-1},dst_comp) = dst({ix,0}, dst_comp);
			dst({this->Nx-1,0},dst_comp) = dst({0,0}, dst_comp);
			#pragma omp parallel for
			for(int ix=0; ix<this->Nx; ix++)
				dst({ix,this->Ny-1},dst_comp) = dst({ix,0}, dst_comp);
		}
		else
			dst.set_block_to_zero(dst_comp);
	}
}

template <typename T, int S>
void BlockSparseMatrixDXDY<T,S>::shifted_block_vmult(
		std::complex<double> scaling_factor, const Block &block_idx,
		TransverseOpticalField &dst, int dst_comp,
		TransverseOpticalField &src, int src_comp,
		bool add) const {

	std::complex<double> val;
	if(add) {
		if(this->Nx>3 && this->Ny>3) {
			dst().segment(dst_comp*this->N,this->N) +=
				src().segment(src_comp*this->N,this->N) +
				scaling_factor*this->matrix.block(
					(block_idx/2)*this->N,(block_idx%2)*this->N,this->N,this->N)*
						src().segment(src_comp*this->N,this->N);
			#pragma omp parallel for private(val)
			for(int iy=1; iy<this->Ny-1; iy++) {
				// int ix = 0;
				// val =
					// (*this)(block_idx,{ix,iy},pxpyShift) * src({ix+1,iy+1},src_comp) +
					// (*this)(block_idx,{ix,iy},pxmyShift) * src({ix+1,iy-1},src_comp) +
					// (*this)(block_idx,{ix,iy},mxpyShift) * src({this->Nx-2,iy+1},src_comp) +
					// (*this)(block_idx,{ix,iy},mxmyShift) * src({this->Nx-2,iy-1},src_comp);
				// dst({ix,iy},dst_comp) +=
					// src({ix,iy},src_comp) + scaling_factor*val;
				// for(ix=1; ix<this->Nx-1; ix++) {
					// val =
						// (*this)(block_idx,{ix,iy},pxpyShift) * src({ix+1,iy+1},src_comp) +
						// (*this)(block_idx,{ix,iy},pxmyShift) * src({ix+1,iy-1},src_comp) +
						// (*this)(block_idx,{ix,iy},mxpyShift) * src({ix-1,iy+1},src_comp) +
						// (*this)(block_idx,{ix,iy},mxmyShift) * src({ix-1,iy-1},src_comp);
					// dst({ix,iy},dst_comp) +=
						// src({ix,iy},src_comp) + scaling_factor*val;
				// }
				// dst({ix,iy},dst_comp) = dst({0,iy}, dst_comp);
				dst({this->Nx-1,iy},dst_comp) = dst({0,iy}, dst_comp);
			}
			// int ix = 0;
			// val =
				// (*this)(block_idx,{ix,0},pxpyShift) * src({ix+1,1},src_comp) +
				// (*this)(block_idx,{ix,0},pxmyShift) * src({ix+1,this->Ny-2},src_comp) +
				// (*this)(block_idx,{ix,0},mxpyShift) * src({this->Nx-2,1},src_comp) +
				// (*this)(block_idx,{ix,0},mxmyShift) * src({this->Nx-2,this->Ny-2},src_comp);
			// dst({ix,0},dst_comp) +=
				// src({ix,0},src_comp) + scaling_factor*val;
			// for(ix=1; ix<this->Nx-1; ix++) {
				// val =
					// (*this)(block_idx,{ix,0},pxpyShift) * src({ix+1,1},src_comp) +
					// (*this)(block_idx,{ix,0},pxmyShift) * src({ix+1,+this->Ny-2},src_comp) +
					// (*this)(block_idx,{ix,0},mxpyShift) * src({ix-1,1},src_comp) +
					// (*this)(block_idx,{ix,0},mxmyShift) * src({ix-1,+this->Ny-2},src_comp);
				// dst({ix,0},dst_comp) +=
					// src({ix,0},src_comp) + scaling_factor*val;
			// }
			// dst({ix,0},dst_comp) = dst({0,0}, dst_comp);
			// for(ix=0; ix<this->Nx; ix++)
				// dst({ix,this->Ny-1},dst_comp) = dst({ix,0}, dst_comp);
			dst({this->Nx-1,0},dst_comp) = dst({0,0}, dst_comp);
			#pragma omp parallel for
			for(int ix=0; ix<this->Nx; ix++)
				dst({ix,this->Ny-1},dst_comp) = dst({ix,0}, dst_comp);
		}
		else
			dst.add_block(dst_comp, src, src_comp);
	}
	else {
		if(this->Nx>3 && this->Ny>3) {
			dst().segment(dst_comp*this->N,this->N) =
				src().segment(src_comp*this->N,this->N) +
				scaling_factor*this->matrix.block(
					(block_idx/2)*this->N,(block_idx%2)*this->N,this->N,this->N)*
						src().segment(src_comp*this->N,this->N);
			#pragma omp parallel for private(val)
			for(int iy=1; iy<this->Ny-1; iy++) {
				// int ix = 0;
				// val =
					// (*this)(block_idx,{ix,iy},pxpyShift) * src({ix+1,iy+1},src_comp) +
					// (*this)(block_idx,{ix,iy},pxmyShift) * src({ix+1,iy-1},src_comp) +
					// (*this)(block_idx,{ix,iy},mxpyShift) * src({this->Nx-2,iy+1},src_comp) +
					// (*this)(block_idx,{ix,iy},mxmyShift) * src({this->Nx-2,iy-1},src_comp);
				// dst({ix,iy},dst_comp) =
					// src({ix,iy},src_comp) + scaling_factor*val;
				// for(ix=1; ix<this->Nx-1; ix++) {
					// val =
						// (*this)(block_idx,{ix,iy},pxpyShift) * src({ix+1,iy+1},src_comp) +
						// (*this)(block_idx,{ix,iy},pxmyShift) * src({ix+1,iy-1},src_comp) +
						// (*this)(block_idx,{ix,iy},mxpyShift) * src({ix-1,iy+1},src_comp) +
						// (*this)(block_idx,{ix,iy},mxmyShift) * src({ix-1,iy-1},src_comp);
					// dst({ix,iy},dst_comp) =
						// src({ix,iy},src_comp) + scaling_factor*val;
				// }
				// dst({ix,iy},dst_comp) = dst({0,iy}, dst_comp);
				dst({this->Nx-1,iy},dst_comp) = dst({0,iy}, dst_comp);
			}
			// int ix = 0;
			// val =
				// (*this)(block_idx,{ix,0},pxpyShift) * src({ix+1,1},src_comp) +
				// (*this)(block_idx,{ix,0},pxmyShift) * src({ix+1,this->Ny-2},src_comp) +
				// (*this)(block_idx,{ix,0},mxpyShift) * src({this->Nx-2,1},src_comp) +
				// (*this)(block_idx,{ix,0},mxmyShift) * src({this->Nx-2,this->Ny-2},src_comp);
			// dst({ix,0},dst_comp) =
				// src({ix,0},src_comp) + scaling_factor*val;
			// for(ix=1; ix<this->Nx-1; ix++) {
				// val =
					// (*this)(block_idx,{ix,0},pxpyShift) * src({ix+1,1},src_comp) +
					// (*this)(block_idx,{ix,0},pxmyShift) * src({ix+1,+this->Ny-2},src_comp) +
					// (*this)(block_idx,{ix,0},mxpyShift) * src({ix-1,1},src_comp) +
					// (*this)(block_idx,{ix,0},mxmyShift) * src({ix-1,+this->Ny-2},src_comp);
				// dst({ix,0},dst_comp) =
					// src({ix,0},src_comp) + scaling_factor*val;
			// }
			// dst({ix,0},dst_comp) = dst({0,0}, dst_comp);
			// for(ix=0; ix<this->Nx; ix++)
				// dst({ix,this->Ny-1},dst_comp) = dst({ix,0}, dst_comp);
			dst({this->Nx-1,0},dst_comp) = dst({0,0}, dst_comp);
			#pragma omp parallel for
			for(int ix=0; ix<this->Nx; ix++)
				dst({ix,this->Ny-1},dst_comp) = dst({ix,0}, dst_comp);
		}
		else
			dst.copy_block(dst_comp, src, src_comp);
	}
}

template <typename T, int S>
SplittedBlockSparseMatrix<T,S>::SplittedBlockSparseMatrix(
		const CartesianMesh &mesh,
		const bool (&zero_blocks_DX)[4],
		const bool (&zero_blocks_DY)[4],
		const bool (&zero_blocks_DXDY)[4]) :
	DX_mat(mesh, zero_blocks_DX),
	DY_mat(mesh, zero_blocks_DY),
	DXDY_mat(mesh, zero_blocks_DXDY),
	mats{&DX_mat, &DY_mat, &DXDY_mat} {}

template <typename T, int S>
void SplittedBlockSparseMatrix<T,S>::block_vmult(
		const Block &block_idx,
		TransverseOpticalField &dst, int dst_comp,
		TransverseOpticalField &src, int src_comp,
		bool add) const {

	DX_mat.block_vmult(block_idx, dst, dst_comp, src, src_comp, add);
	DY_mat.block_vmult(block_idx, dst, dst_comp, src, src_comp, true);
	DXDY_mat.block_vmult(block_idx, dst, dst_comp, src, src_comp, true);
}

template class BlockSparseMatrix<double,Eigen::ColMajor>;
template class BlockSparseMatrix<double,Eigen::RowMajor>;
template class BlockSparseMatrix<std::complex<double>,Eigen::ColMajor>;
template class BlockSparseMatrix<std::complex<double>,Eigen::RowMajor>;

template class BlockSparseMatrixDX<double,Eigen::ColMajor>;
template class BlockSparseMatrixDX<double,Eigen::RowMajor>;
template class BlockSparseMatrixDX<std::complex<double>,Eigen::ColMajor>;
template class BlockSparseMatrixDX<std::complex<double>,Eigen::RowMajor>;

template class BlockSparseMatrixDY<double,Eigen::ColMajor>;
template class BlockSparseMatrixDY<double,Eigen::RowMajor>;
template class BlockSparseMatrixDY<std::complex<double>,Eigen::ColMajor>;
template class BlockSparseMatrixDY<std::complex<double>,Eigen::RowMajor>;

template class BlockSparseMatrixDXDY<double,Eigen::ColMajor>;
template class BlockSparseMatrixDXDY<double,Eigen::RowMajor>;
template class BlockSparseMatrixDXDY<std::complex<double>,Eigen::ColMajor>;
template class BlockSparseMatrixDXDY<std::complex<double>,Eigen::RowMajor>;

template class SplittedBlockSparseMatrix<double,Eigen::ColMajor>;
template class SplittedBlockSparseMatrix<double,Eigen::RowMajor>;
template class SplittedBlockSparseMatrix<std::complex<double>,Eigen::ColMajor>;
template class SplittedBlockSparseMatrix<std::complex<double>,Eigen::RowMajor>;
