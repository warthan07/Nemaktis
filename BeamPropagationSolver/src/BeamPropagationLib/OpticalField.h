#ifndef OPTICALFIELD_H
#define OPTICALFIELD_H

#include <Eigen/SparseCore>
#include "PermittivityTensorField.h"

class TransverseOpticalField {
public:
	TransverseOpticalField(const CartesianMesh &mesh) :
		Nx(mesh.Nx),
		Ny(mesh.Ny),
		N(Nx*Ny),
		vec(2*Nx*Ny) {}

	void operator=(const TransverseOpticalField &src) {
		vec = src.vec;
	}
	void operator=(const Eigen::VectorXcd &src) {
		vec = src;
	}

	void add_scaled_field(std::complex<double> &scaling_factor, TransverseOpticalField &src) {
		vec += scaling_factor*src.vec;
	}
	inline void set_block_to_zero(int comp) {
		vec.segment(comp*N, N).setZero();
	}
	inline void add_block(
			const TransverseOpticalField &src, int comp) {
		vec.segment(comp*N, N) += src().segment(comp*N, N);
	}
	inline void add_block(
			int dst_comp, const TransverseOpticalField &src, int src_comp) {
		vec.segment(dst_comp*N, N) += src().segment(src_comp*N, N);
	}
	inline void copy_block(
			const TransverseOpticalField &src, int comp) {
		vec.segment(comp*N, N) = src().segment(comp*N, N);
	}
	inline void copy_block(
			int dst_comp, const TransverseOpticalField &src, int src_comp) {
		vec.segment(dst_comp*N, N) = src().segment(src_comp*N, N);
	}
	inline void scale_block(
			const std::complex<double> &scaling_factor, int comp) {
		vec.segment(comp*N, N) *= scaling_factor;
	}
	inline void scale_and_add_block(
			const std::complex<double> &scaling_factor,
			TransverseOpticalField &src, int comp) {
		vec.segment(comp*N, N) *= scaling_factor;
		vec.segment(comp*N, N) += src().segment(comp*N, N);
	}

	const Eigen::VectorXcd& operator()() const {
		return vec;
	}
	Eigen::VectorXcd& operator()() {
		return vec;
	}

	std::complex<double> operator()(int mesh_idx, int comp) const {
		return vec[mesh_idx+N*comp];
	}
	std::complex<double>& operator()(int mesh_idx, int comp) {
		return vec[mesh_idx+N*comp];
	}

	std::complex<double> operator()(const Index2D &p, int comp) const {
		return vec[p.x+Nx*p.y+N*comp];
	}
	std::complex<double>& operator()(const Index2D &p, int comp) {
		return vec[p.x+Nx*p.y+N*comp];
	}

private:
	const int Nx, Ny, N;
	Eigen::VectorXcd vec; 
};

class LongitudinalOpticalField {
public:
	LongitudinalOpticalField(const CartesianMesh &mesh) :
		Nx(mesh.Nx),
		Ny(mesh.Ny),
		vec(2*Nx*Ny) {}

	void operator=(const LongitudinalOpticalField &src) {
		vec = src.vec;
	}
	void operator=(const Eigen::VectorXcd &src) {
		vec = src;
	}

	const Eigen::VectorXcd& operator()() const {
		return vec;
	}
	Eigen::VectorXcd& operator()() {
		return vec;
	}

	std::complex<double> operator()(int mesh_idx) const {
		return vec[mesh_idx];
	}
	std::complex<double>& operator()(int mesh_idx) {
		return vec[mesh_idx];
	}

	std::complex<double> operator()(const Index2D &p) const {
		return vec[p.x+Nx*p.y];
	}
	std::complex<double>& operator()(const Index2D &p) {
		return vec[p.x+Nx*p.y];
	}

private:
	const int Nx, Ny;
	Eigen::VectorXcd vec; 
};


#endif
