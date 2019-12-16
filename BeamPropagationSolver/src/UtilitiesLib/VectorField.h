#ifndef VECTORFIELD_H
#define VECTORFIELD_H

#include <vector>
#include <string>

#include "CartesianMesh.h"

struct Index2D {
	/**
	 * Index for each spatial directions.
	 */
	int x, y;
};


struct Index3D {
	/**
	 * Index for each spatial directions.
	 */
	int x, y, z;
};

template <typename T>
class VectorField {
public:
	VectorField(CartesianMesh mesh, int field_dim);
	VectorField(CartesianMesh mesh, int field_dim, T val);
	/**
	 * Initializes a vector field from a raw pointer array (slow to fast
	 * indices: iz -> iy -> ix -> field component). The number of components
	 * should either be field_dim (director field without any
	 * mask) or field_dim+1 (director field with a mask specified in the
	 * last component). Positive mask values are equivalent to true,
	 * strictly negative mask values are equivalent to false.
	 */
	VectorField(CartesianMesh mesh, int field_dim, T* user_vals, int n_user_vals);

	void set_mask(std::string mask_formula);
	void set_mask_from_nonzeros();
	bool get_mask_val(const int vertex_idx) const;
	bool get_mask_val(const Index3D &p) const;

	void operator=(T val);
	void operator=(const VectorField<T> &src);

	void operator*=(T val);
	void operator/=(T val);

	void operator+=(const VectorField<T> &src);
	void operator-=(const VectorField<T> &src);

	T operator[](int idx) const;
	T& operator[](int idx);

	T operator()(const int vertex_idx, int comp) const;
	T& operator()(const int vertex_idx, int comp);

	T operator()(const Index3D &p, int comp) const;
	T& operator()(const Index3D &p, int comp);

	T operator,(const VectorField<T> &src);

	void add(T scaling_factor, const VectorField<T> &src);

	double norm_sqr();
	double norm();
	double linfty_norm();

	const int field_dim;
	const int n_vertex;

	const CartesianMesh mesh;

protected:
	std::vector<T> vals;
	std::vector<bool> mask_vals;

	bool mask_exists;
};

#endif
