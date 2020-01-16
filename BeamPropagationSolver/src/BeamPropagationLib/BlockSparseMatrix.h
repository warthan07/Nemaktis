#ifndef BLOCKSPARSEMATRIX_H
#define BLOCKSPARSEMATRIX_H

#include <memory>
#include <Eigen/SparseCore>

#include "CartesianMesh.h"
#include "PermittivityTensorField.h"
#include "OpticalField.h"

/**
 * Enum for acessing the sub-blocks of a 2x2 block matrix
 */
enum Block {
	Block00 = 0,
	Block01 = 1,
	Block10 = 2,
	Block11 = 3
};

/**
 * Enum for shifting a mesh index to a nearby mesh point
 */
enum Shift {
	noShift = 0,
	pxShift = 1,
	mxShift = 2,
	pyShift = 3,
	myShift = 4,
	pxpyShift = 5,
	mxpyShift = 6,
	pxmyShift = 7,
	mxmyShift = 8
};

/**
 * Enum for selecting the type of differential operator: pure
 * x-derivatives, pure y-derivatives, mixed xy-derivatives 
 */
enum DiffType {
	DiffDX,
	DiffDY,
	DiffDXDY
};

/**
 * A wrapper class of Eigen::SparseMatrix allowing quick O(1) access to
 * matrix elements based on a mesh index, block index, and "shift"
 * vector on the stencil, for a 2D differential operator (stencil is
 * therefore included on 3x3 2D grid). With this base object, it is the
 * user responsability to initialize the sparsity pattern of the
 * underlying matrix, and then call the method update_pointers() before
 * using the access method for the matrix entry.
 */
template <typename T, int Scheme=Eigen::ColMajor>
class BlockSparseMatrix {
public:
	BlockSparseMatrix(CartesianMesh mesh);

	void update_pointers();

	Eigen::SparseMatrix<T,Scheme>& operator()() {
		return matrix;
	}
	const Eigen::SparseMatrix<T,Scheme>& operator()() const {
		return matrix;
	}
	T& operator()(
			const Block &block_idx, const Index2D &mesh_idx, const Shift &shift_idx) {
		return *ptrs[block_idx][shift_idx][mesh_idx.x+Nx*mesh_idx.y];
	}
	T operator()(
			const Block &block_idx, const Index2D &mesh_idx, const Shift &shift_idx) const {
		return *ptrs[block_idx][shift_idx][mesh_idx.x+Nx*mesh_idx.y];
	}

	struct ColIteratorEntry {
		Eigen::Index row;
		T& val;
	};
	struct RowIteratorEntry {
		Eigen::Index col;
		T& val;
	};

	std::vector<ColIteratorEntry>& column_iterators(int col_index) {
		return col_iterator_arrays[col_index];
	}
	std::vector<RowIteratorEntry>& row_iterators(int row_index) {
		return row_iterator_arrays[row_index];
	}

protected:
	const int Nx, Ny, N;
	Eigen::SparseMatrix<T,Scheme> matrix;
	std::vector<T*> ptrs[4][9];
	std::vector<std::vector<ColIteratorEntry> > col_iterator_arrays;
	std::vector<std::vector<RowIteratorEntry> > row_iterator_arrays;
};

/**
 * A more specialized version BlockSparseMatrix where each sub-block is
 * a tridiagonal matrix corresponding to derivatives of the type d/dx
 * and d/dx^2. Additional methods are given for sub-block matrix-vector
 * multiplication and inverse. Here, the sparsity pattern is
 * automatically built at construction and the update_pointers method of
 * the parent class is automatically called.
 */
template <typename T, int Scheme=Eigen::ColMajor>
class BlockSparseMatrixDX : public BlockSparseMatrix<T,Scheme> {
public:
	BlockSparseMatrixDX(
		const CartesianMesh &mesh, 
		const bool (&zero_blocks)[4] = {false, false, false, false});

	/**
	 * For a given sub-block A=mat[i][j] of the matrix of this object,
	 * apply the following matrix-vector product operation:
	 * dst[dst_comp] = A * src[src_comp] if add=false,
	 * dst[dst_comp] += A * src[src_comp] if add=true,
	 *
	 */
	void block_vmult(
			const Block &block_idx,
			TransverseOpticalField &dst, int dst_comp,
			TransverseOpticalField &src, int src_comp,
			bool add = false) const;
	/**
	 * For a given sub-block A=mat[i][j] of the matrix of this object,
	 * apply the following matrix-vector product operation:
	 * dst[dst_comp] = (Id+scaling_factor*A) * src[src_comp] if add=false,
	 * dst[dst_comp] += (Id+scaling_factor*A) * src[src_comp] if add=true,
	 */
	void shifted_block_vmult(
			std::complex<double> scaling_factor, const Block &block_idx,
			TransverseOpticalField &dst, int dst_comp,
			TransverseOpticalField &src, int src_comp,
			bool add = false) const;
	/**
	 * For a given sub-block A=mat[i][j] of the matrix of this object,
	 * apply the following inverse matrix-vector product operation:
	 * dst[dst_comp] = (Id+scaling_factor*A)^-1 * src[src_comp].
	 */
	void shifted_block_vmult_inv(
			std::complex<double> scaling_factor, const Block &block_idx,
			TransverseOpticalField &dst, int dst_comp,
			const TransverseOpticalField &src, int src_comp) const;

private:
	inline std::complex<double> diag_sup(
			std::complex<double> scaling_coef, const Block &block_idx, const Index2D &p) const {
		return scaling_coef*(*this)(block_idx, p, pxShift);
	}
	inline std::complex<double> diag_mid(
			std::complex<double> scaling_coef, const Block &block_idx, const Index2D &p) const {
		return 1.+scaling_coef*(*this)(block_idx, p, noShift);
	}
	inline std::complex<double> diag_inf(
			std::complex<double> scaling_coef, const Block &block_idx, const Index2D &p) const {
		return scaling_coef*(*this)(block_idx, p, mxShift);
	}
};

/**
 * A more specialized version BlockSparseMatrix where each sub-block is
 * a tridiagonal matrix corresponding to derivatives of the type d/dy
 * and d/dy^2. Additional methods are given for sub-block matrix-vector
 * multiplication and inverse. Here, the sparsity pattern is
 * automatically built at construction and the update_pointers method of
 * the parent class is automatically called.
 */
template <typename T, int Scheme=Eigen::ColMajor>
class BlockSparseMatrixDY : public BlockSparseMatrix<T,Scheme> {
public:
	BlockSparseMatrixDY(
		const CartesianMesh &mesh,
		const bool (&zero_blocks)[4] = {false, false, false, false});

	/**
	 * For a given sub-block A=mat[i][j] of the matrix of this object,
	 * apply the following matrix-vector product operation:
	 * dst[dst_comp] = A * src[src_comp] if add=false,
	 * dst[dst_comp] += A * src[src_comp] if add=true,
	 *
	 */
	void block_vmult(
			const Block &block_idx,
			TransverseOpticalField &dst, int dst_comp,
			TransverseOpticalField &src, int src_comp,
			bool add = false) const;
	/**
	 * For a given sub-block A=mat[i][j] of the matrix of this object,
	 * apply the following matrix-vector product operation:
	 * dst[dst_comp] = (Id+scaling_factor*A) * src[src_comp] if add=false,
	 * dst[dst_comp] += (Id+scaling_factor*A) * src[src_comp] if add=true,
	 */
	void shifted_block_vmult(
			std::complex<double> scaling_factor, const Block &block_idx,
			TransverseOpticalField &dst, int dst_comp,
			TransverseOpticalField &src, int src_comp,
			bool add = false) const;
	/**
	 * For a given sub-block A=mat[i][j] of the matrix of this object,
	 * apply the following inverse matrix-vector product operation:
	 * dst[dst_comp] = (Id+scaling_factor*A)^-1 * src[src_comp].
	 */
	void shifted_block_vmult_inv(
			std::complex<double> scaling_factor, const Block &block_idx,
			TransverseOpticalField &dst, int dst_comp,
			const TransverseOpticalField &src, int src_comp) const;

private:
	inline std::complex<double> diag_sup(
			std::complex<double> scaling_coef, const Block &block_idx, const Index2D &p) const {
		return scaling_coef*(*this)(block_idx, p, pyShift);
	}
	inline std::complex<double> diag_mid(
			std::complex<double> scaling_coef, const Block &block_idx, const Index2D &p) const {
		return 1.+scaling_coef*(*this)(block_idx, p, noShift);
	}
	inline std::complex<double> diag_inf(
			std::complex<double> scaling_coef, const Block &block_idx, const Index2D &p) const {
		return scaling_coef*(*this)(block_idx, p, myShift);
	}
};

/**
 * A more specialized version BlockSparseMatrix where each sub-block is
 * a quadridiagonal matrix corresponding to derivatives of the type
 * d/dxdy. Additional methods are given for sub-block matrix-vector
 * multiplication. Here, the sparsity pattern is automatically built at
 * construction and the update_pointers method of the parent class is
 * automatically called.
 */
template <typename T, int Scheme=Eigen::ColMajor>
class BlockSparseMatrixDXDY : public BlockSparseMatrix<T,Scheme> {
public:
	BlockSparseMatrixDXDY(
		const CartesianMesh &mesh,
		const bool (&zero_blocks)[4] = {false, false, false, false});

	/**
	 * For a given sub-block A=mat[i][j] of the matrix of this object,
	 * apply the following matrix-vector product operation:
	 * dst[dst_comp] = A * src[src_comp] if add=false,
	 * dst[dst_comp] += A * src[src_comp] if add=true,
	 *
	 */
	void block_vmult(
			const Block &block_idx,
			TransverseOpticalField &dst, int dst_comp,
			TransverseOpticalField &src, int src_comp,
			bool add = false) const;
	/**
	 * For a given sub-block A=mat[i][j] of the matrix of this object,
	 * apply the following matrix-vector product operation:
	 * dst[dst_comp] = (Id+scaling_factor*A) * src[src_comp] if add=false,
	 * dst[dst_comp] += (Id+scaling_factor*A) * src[src_comp] if add=true,
	 */
	void shifted_block_vmult(
			std::complex<double> scaling_factor, const Block &block_idx,
			TransverseOpticalField &dst, int dst_comp,
			TransverseOpticalField &src, int src_comp,
			bool add = false) const;
};

template <typename T, int S=Eigen::ColMajor>
class SplittedBlockSparseMatrix {
public:
	SplittedBlockSparseMatrix(
		const CartesianMesh &mesh,
		const bool (&zero_blocks_DX)[4] = {false, false, false, false},
		const bool (&zero_blocks_DY)[4] = {false, false, false, false},
		const bool (&zero_blocks_DXDY)[4] = {false, false, false, false});

	void block_vmult(
			const Block &block_idx,
			TransverseOpticalField &dst, int dst_comp,
			TransverseOpticalField &src, int src_comp,
			bool add = false) const;

	T& operator()(
			const DiffType &diff_type, const Block &block_idx,
			const Index2D &mesh_idx, const Shift &shift_idx) {
		return (*mats[diff_type])(block_idx, mesh_idx, shift_idx);
	}
	T operator()(
			const DiffType &diff_type, const Block &block_idx,
			const Index2D &mesh_idx, const Shift &shift_idx) const {
		return (*mats[diff_type])(block_idx, mesh_idx, shift_idx);
	}

	std::vector<typename BlockSparseMatrix<T,S>::ColIteratorEntry>& column_iterators(
			const DiffType &diff_type, int col_index) {
		return mats[diff_type]->column_iterators(col_index);
	}
	std::vector<typename BlockSparseMatrix<T,S>::RowIteratorEntry>& row_iterators(
			const DiffType &diff_type, int row_index) {
		return mats[diff_type]->row_iterators(row_index);
	}

	const BlockSparseMatrixDX<T,S>& DX_op() {
		return DX_mat;
	}
	const BlockSparseMatrixDY<T,S>& DY_op() {
		return DY_mat;
	}
	const BlockSparseMatrixDXDY<T,S>& DXDY_op() {
		return DXDY_mat;
	}

private:
	BlockSparseMatrixDX<T,S> DX_mat;
	BlockSparseMatrixDY<T,S> DY_mat;
	BlockSparseMatrixDXDY<T,S> DXDY_mat;

	BlockSparseMatrix<T,S>* mats[3];
};


#endif
