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
	inline void add_block(
			TransverseOpticalField &src, int comp) {
		#pragma omp parallel for
		for(int iy=1; iy<Ny-1; iy++) 
			for(int ix=1; ix<Nx-1; ix++) 
				(*this)({ix,iy},comp) += src({ix,iy},comp);
	}
	inline void copy_block(
			TransverseOpticalField &src, int comp) {
		#pragma omp parallel for
		for(int iy=1; iy<Ny-1; iy++) 
			for(int ix=1; ix<Nx-1; ix++) 
				(*this)({ix,iy},comp) = src({ix,iy},comp);
	}
	inline void scale_block(
			const std::complex<double> &scaling_factor, int comp) {
		#pragma omp parallel for
		for(int iy=1; iy<Ny-1; iy++) 
			for(int ix=1; ix<Nx-1; ix++) 
				(*this)({ix,iy},comp) *= scaling_factor;
	}
	inline void scale_and_add_block(
			const std::complex<double> &scaling_factor,
			TransverseOpticalField &src, int comp) {
		#pragma omp parallel for
		for(int iy=1; iy<Ny-1; iy++) 
			for(int ix=1; ix<Nx-1; ix++) 
				(*this)({ix,iy},comp) = src({ix,iy},comp) + scaling_factor*(*this)({ix,iy},comp);
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
