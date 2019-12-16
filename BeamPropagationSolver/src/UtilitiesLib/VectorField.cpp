#include <complex>
#include <algorithm>
#include <muParser.h>

#include "VectorField.h"
#include "error.h"

template <typename T>
VectorField<T>::VectorField(CartesianMesh mesh, int field_dim) :
	field_dim(field_dim),
	n_vertex(mesh.Nx*mesh.Ny*mesh.Nz),
	mesh(mesh),
	vals(n_vertex * field_dim),
	mask_exists(false) {}

template <typename T>
VectorField<T>::VectorField(CartesianMesh mesh, int field_dim, T val) :
	field_dim(field_dim),
	n_vertex(mesh.Nx*mesh.Ny*mesh.Nz),
	mesh(mesh),
	vals(n_vertex * field_dim, val),
	mask_exists(false) {}

template <typename T>
VectorField<T>::VectorField(CartesianMesh mesh, int field_dim, T* user_vals, int n_user_vals) :
	field_dim(field_dim),
	n_vertex(mesh.Nx*mesh.Ny*mesh.Nz),
	mesh(mesh),
	vals(n_vertex * field_dim),
	mask_exists(false) {

	if(n_user_vals!=n_vertex*field_dim && n_user_vals!=n_vertex*(field_dim+1))
		throw std::string("Wrong dimension for the raw pointer array");

	if(n_user_vals==n_vertex*field_dim)
		for(int i=0; i<n_vertex*field_dim; i++)
			vals[i] = user_vals[i];
	else {
		mask_exists = true;
		mask_vals.resize(n_vertex);
		for(int i=0; i<n_vertex; i++) {
			for(int c=0; c<field_dim; c++)
				(*this)(i,c) = user_vals[(field_dim+1)*i+c];
			mask_vals[i] = (std::real(user_vals[(field_dim+1)*i+field_dim])>=0);
		}
	}
}

template <typename T>
void VectorField<T>::set_mask(std::string mask_formula) {

	double x, y, z;
	mask_vals.resize(n_vertex);
	mask_exists = true;

	#pragma omp parallel for firstprivate(x, y, z)
	for(int iz=0; iz<mesh.Nz; iz++) {
		z = (iz-0.5*(mesh.Nz-1))*mesh.delta_z;

		mu::Parser p;
		p.DefineVar("x", &x);
		p.DefineVar("y", &y);
		p.DefineConst("z", z);
		p.SetExpr(mask_formula);

		for(int iy=0; iy<mesh.Ny; iy++) {
			y = (iy-0.5*(mesh.Ny-1))*mesh.delta_y;
			for(int ix=0; ix<mesh.Nx; ix++) {
				x = (ix-0.5*(mesh.Nx-1))*mesh.delta_x;
				mask_vals[ix + mesh.Nx*(iy + mesh.Ny*iz)] = (p.Eval()>=0);
			}
		}
	}
}

template <typename T>
void VectorField<T>::set_mask_from_nonzeros() {

	mask_vals.resize(n_vertex);
	mask_exists = true;
	#pragma omp parallel for
	for(int vertex_idx=0; vertex_idx<n_vertex; vertex_idx++) {
		double norm = 0;
		for(int c=0; c<field_dim; c++)
			norm += std::pow(std::abs((*this)(vertex_idx,c)), 2);
		norm = std::sqrt(norm);
		mask_vals[vertex_idx] = (norm>1e-8);
	}
}

template <typename T>
bool VectorField<T>::get_mask_val(const int vertex_idx) const {

	if(!mask_exists)
		return true;
	else
		return mask_vals[vertex_idx];
}

template <typename T>
bool VectorField<T>::get_mask_val(const Index3D &p) const {

	if(!mask_exists)
		return true;
	else
		return mask_vals[p.x + mesh.Nx*( p.y + mesh.Ny*p.z)];
}

template <typename T>
void VectorField<T>::operator=(T val) {
	
	for(int i=0; i<vals.size(); i++)
		vals[i] = val;
}

template <typename T>
void VectorField<T>::operator=(const VectorField<T> &src) {

	Assert(
		src.n_vertex == n_vertex && src.field_dim == field_dim,
		"Cannot copy vector field: wrong dimension");
	
	for(int i=0; i<vals.size(); i++)
		vals[i] = src[i];
}

template <typename T>
void VectorField<T>::operator*=(T val) {

	#pragma omp parallel for
	for(int iz=0; iz<mesh.Nz; iz++)
		for(int iy=0; iy<mesh.Ny; iy++)
			for(int ix=0; ix<mesh.Nx; ix++)
				for(int comp=0; comp<field_dim; comp++)
					(*this)({ix,iy,iz}, comp) *= val;
}

template <typename T>
void VectorField<T>::operator/=(T val) {

	#pragma omp parallel for
	for(int iz=0; iz<mesh.Nz; iz++)
		for(int iy=0; iy<mesh.Ny; iy++)
			for(int ix=0; ix<mesh.Nx; ix++)
				for(int comp=0; comp<field_dim; comp++)
					(*this)({ix,iy,iz}, comp) /= val;
}

template <typename T>
void VectorField<T>::operator+=(const VectorField<T> &src) {

	#pragma omp parallel for
	for(int iz=0; iz<mesh.Nz; iz++)
		for(int iy=0; iy<mesh.Ny; iy++)
			for(int ix=0; ix<mesh.Nx; ix++)
				for(int comp=0; comp<field_dim; comp++)
					(*this)({ix,iy,iz}, comp) += src({ix,iy,iz}, comp);
}

template <typename T>
void VectorField<T>::operator-=(const VectorField<T> &src) {

	#pragma omp parallel for
	for(int iz=0; iz<mesh.Nz; iz++)
		for(int iy=0; iy<mesh.Ny; iy++)
			for(int ix=0; ix<mesh.Nx; ix++)
				for(int comp=0; comp<field_dim; comp++)
					(*this)({ix,iy,iz}, comp) -= src({ix,iy,iz}, comp);
}

template <typename T>
T VectorField<T>::operator[](int idx) const {

	Assert(
		idx<n_vertex * field_dim, "Out-of-range index (VectorField)");

	return vals[idx];
}

template <typename T>
T& VectorField<T>::operator[](int idx) {

	Assert(
		idx<n_vertex * field_dim, "Out-of-range index (VectorField)");

	return vals[idx];
}

template <typename T>
T VectorField<T>::operator()(const int vertex_idx, int comp) const {

	Assert(
		vertex_idx<n_vertex && comp<field_dim,
		"Out of range index (VectorField)");

	return vals[field_dim*vertex_idx + comp];
}

template <typename T>
T& VectorField<T>::operator()(const int vertex_idx, int comp) {

	Assert(
		vertex_idx<n_vertex && comp<field_dim,
		"Out of range index (VectorField)");

	return vals[field_dim*vertex_idx + comp];
}

template <typename T>
T VectorField<T>::operator()(const Index3D &p, int comp) const {

	Assert(
		p.x<mesh.Nx && p.y<mesh.Ny && p.z<mesh.Nz && comp<field_dim,
		"Out-of-range index (VectorField)");

	return vals[comp + field_dim*( p.x + mesh.Nx*( p.y + mesh.Ny*p.z ) )];
}

template <typename T>
T& VectorField<T>::operator()(const Index3D &p, int comp) {

	Assert(
		p.x<mesh.Nx && p.y<mesh.Ny && p.z<mesh.Nz && comp<field_dim,
		"Out-of-range index (VectorField)");

	return vals[comp + field_dim*( p.x + mesh.Nx*( p.y + mesh.Ny*p.z ) )];
}

template <>
double VectorField<double>::operator,(const VectorField<double> &src) {

	double res = 0;

	#pragma omp parallel for reduction(+:res)
	for(int iz=0; iz<mesh.Nz; iz++)
		for(int iy=0; iy<mesh.Ny; iy++)
			for(int ix=0; ix<mesh.Nx; ix++)
				for(int comp=0; comp<field_dim; comp++)
					res += (*this)({ix,iy,iz}, comp) * src({ix,iy,iz}, comp);

	return res;
}

template <>
std::complex<double> VectorField<std::complex<double> >::operator,(
		const VectorField<std::complex<double> > &src) {

	double res_re = 0;
	double res_im = 0;

	#pragma omp parallel for reduction(+:res_re,res_im)
	for(int iz=0; iz<mesh.Nz; iz++) {
		for(int iy=0; iy<mesh.Ny; iy++) {
			for(int ix=0; ix<mesh.Nx; ix++) {
				for(int comp=0; comp<field_dim; comp++) {
					res_re += std::real((*this)({ix,iy,iz}, comp) * src({ix,iy,iz}, comp));
					res_im += std::imag((*this)({ix,iy,iz}, comp) * src({ix,iy,iz}, comp));
				}
			}
		}
	}

	return std::complex<double>(res_re,res_im);
}

template <typename T>
void VectorField<T>::add(T scaling_factor, const VectorField<T> &src) {

	#pragma omp parallel for
	for(int iz=0; iz<mesh.Nz; iz++)
		for(int iy=0; iy<mesh.Ny; iy++)
			for(int ix=0; ix<mesh.Nx; ix++)
				for(int comp=0; comp<field_dim; comp++)
					(*this)({ix,iy,iz}, comp) += scaling_factor * src({ix,iy,iz}, comp);
}

template <typename T>
double VectorField<T>::norm_sqr() {

	double res = 0;

	#pragma omp parallel for reduction(+:res)
	for(int iz=0; iz<mesh.Nz; iz++)
		for(int iy=0; iy<mesh.Ny; iy++)
			for(int ix=0; ix<mesh.Nx; ix++)
				for(int comp=0; comp<field_dim; comp++)
					res += std::pow(std::abs((*this)({ix,iy,iz}, comp)), 2.);

	return res;
}

template <typename T>
double VectorField<T>::norm() {

	return std::sqrt(norm_sqr());
}

template <typename T>
double VectorField<T>::linfty_norm() {

	double res = 0;

	for(int iz=0; iz<mesh.Nz; iz++)
		for(int iy=0; iy<mesh.Ny; iy++)
			for(int ix=0; ix<mesh.Nx; ix++)
				for(int comp=0; comp<field_dim; comp++)
					res = std::max(res, std::abs((*this)({ix,iy,iz}, comp)));

	return res;
}

template class VectorField<double>;
template class VectorField<std::complex<double> >;
