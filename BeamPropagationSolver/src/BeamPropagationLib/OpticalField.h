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
		own_data(true),
		raw_data(new std::complex<double>[2*Nx*Ny]),
		vec(raw_data, 2*Nx*Ny) {}

	TransverseOpticalField(const CartesianMesh &mesh, std::complex<double>* raw_data) :
		Nx(mesh.Nx),
		Ny(mesh.Ny),
		N(Nx*Ny),
		own_data(false),
		raw_data(raw_data),
		vec(raw_data, 2*Nx*Ny) {}

	~TransverseOpticalField() {
		if(own_data)
			delete raw_data;
	}

	const auto& operator()() const {
		return vec;
	}
	auto& operator()() {
		return vec;
	}

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
	bool own_data;
	std::complex<double>* raw_data;
	Eigen::Map<Eigen::VectorXcd> vec; 
};

#endif
