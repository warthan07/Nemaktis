#ifndef LINALG_H
#define LINALG_H

#define PI 3.141592653589793238

#include <cmath>
#include <complex>
#include <vector>
#include <iostream>
#include <type_traits>

#include "error.h"
#include "json.h"

using json = nlohmann::json;

/**
 * \defgroup linalg Linear Algebra
 * 
 * A collection of classes and function for basic linear algebra
 * operations.
 * 
 * @{
 */

/**
 * Wrapper class for a vector of dimension <tt>dim</tt> and type
 * <tt>T1</tt>.
 */
template <int dim, typename T1>
class Vector {
public:
	/**
	 * Default constructor, which initializes all components to zero.
	 */
	Vector();
	/**
	 * Constant constructor, which initializes all components to <tt>c</tt>.
	 */
	template <typename T2>
	Vector(T2 c);
	/**
	 * Copy constructor with another Vector object
	 */
	template <typename T2>
	Vector(const Vector<dim,T2> &v);
	/**
	 * Copy constructor with a C-Style array.
	 */
	template <typename T2>
	Vector(const T2 (&v)[dim]);
	/**
	 * Copy constructor with an STL vector.
	 */
	template <typename T2>
	Vector(const std::vector<T2> &v);

	/**
	 * Return the l2-norm of this vector.
	 */
	double norm() const;

	/**
	 * Copy operator with another Vector object.
	 */
	template <typename T2>
	Vector<dim,T1>& operator=(const Vector<dim,T2> &v);
	/**
	 * Copy operator with a C-style array.
	 */
	template <typename T2>
	Vector<dim,T1>& operator=(const T2 (&v)[dim]);
	/**
	 * Copy operator with an STL vector.
	 */
	template <typename T2>
	Vector<dim,T1>& operator=(const std::vector<T2> &v);
	/**
	 * Set all components equal to <tt>c</tt>.
	 */
	Vector<dim,T1>& operator=(T1 c);
	/**
	 * Add <tt>v</tt> to this vector.
	 */
	template <typename T2>
	Vector<dim,T1>& operator+=(const Vector<dim,T2> &v);
	/**
	 * Substract <tt>v</tt> to this vector.
	 */
	template <typename T2>
	Vector<dim,T1>& operator-=(const Vector<dim,T2> &v);
	/**
	 * Multiply all components by <tt>c</tt>.
	 */
	template <typename T2>
	Vector<dim,T1>& operator*=(T2 c);
	/**
	 * Divide all components by <tt>c</tt>.
	 */
	template <typename T2>
	Vector<dim,T1>& operator/=(T2 c);

	/**
	 * Indexing operator, supporting read and write operation.
	 */
	T1& operator()(unsigned int i);
	/**
	 * Constant version of the indexing operator, supporting read operations.
	 */
	const T1& operator()(unsigned int i) const;

	template <int dim2, typename T2>
	friend std::ostream& operator<<(std::ostream &os, const Vector<dim2,T2>& v);

protected:
	/**
	 * The low-level data structure storing the vector.
	 */
	T1 data[dim];
};

/**
 * Wrapper class for a matrix of dimension <tt>(dim1,dim2)</tt> and type
 * <tt>T1</tt>.
 */
template <int dim1, int dim2, typename T1>
class Matrix {
public:
	/**
	 * Default constructor, which initializes all components to zero.
	 */
	Matrix();
	/**
	 * Constant constructor, which initializes all components to <tt>c</tt>.
	 */
	Matrix(T1 c);
	/**
	 * Copy constructor with another Matrix object.
	 */
	template <typename T2>
	Matrix(const Matrix<dim1,dim2,T2> &m);
	/**
	 * Copy constructor with a C-Style array.
	 */
	template <typename T2>
	Matrix(const T2 (&m)[dim1*dim2]);

	/**
	 * Return the l2-norm of this matrix.
	 */
	double norm() const;

	/**
	 * Copy operator.
	 */
	template <typename T2>
	Matrix<dim1,dim2,T1>& operator=(const Matrix<dim1,dim2,T2> &m);
	/**
	 * Set all components equal to <tt>c</tt>.
	 */
	Matrix<dim1,dim2,T1>& operator=(T1 c);
	/**
	 * Add <tt>m</tt> to this matrix
	 */
	template <typename T2>
	Matrix<dim1,dim2,T1>& operator+=(const Matrix<dim1,dim2,T2> &m);
	/**
	 * Substract <tt>m</tt> to this matrix
	 */
	template <typename T2>
	Matrix<dim1,dim2,T1>& operator-=(const Matrix<dim1,dim2,T2> &m);
	/**
	 * Multiply all components by <tt>c</tt>.
	 */
	template <typename T2>
	Matrix<dim1,dim2,T1>& operator*=(T2 c);

	/**
	 * Indexing operator, supporting read and write operation.
	 */
	T1& operator()(unsigned int i, unsigned int j);
	/**
	 * Constant version of the indexing operator, supporting read operations.
	 */
	const T1& operator()(unsigned int i, unsigned int j) const;

	template <int dim3, int dim4, typename T2>
	friend std::ostream& operator<<(
		std::ostream &os, const Matrix<dim3,dim4,T2>& m);

protected:
	/**
	 * The low-level data structure storing the vector.
	 */
	T1 data[dim1*dim2];
};

/**
 * A multidimensional index useful for computation and interpolation on
 * n-dimensional meshes.
 */
template <int dim>
class MultiDimIndex {
public:
	MultiDimIndex(
		const Vector<dim,long> &n_points_per_dim,
		const Vector<dim,long> &shift =	Vector<dim,long>((long)0),
		bool last_idx_fastest = true);

	/**
	 * Returns the d-th component of the multidim index.
	 */
	long& operator()(long d);
	const long& operator()(long d) const;
	/**
	 * Returns the flattened index.
	 */
	long operator()() const;
	/**
	 * Returns the flattened index after a shift.
	 */
	long operator()(MultiDimIndex<dim> &shift) const;

	/**
	 * Compute the multidim index associated to a flattened index.
	 */
	MultiDimIndex<dim>& operator=(long idx);
	/**
	 * Copy a multidim index.
	 */
	MultiDimIndex<dim>& operator=(const MultiDimIndex &multi_idx);

	/**
	 * Increment a multidim index.
	 */
	MultiDimIndex<dim>& operator++();
	MultiDimIndex<dim> operator++(int);
	/**
	 * Check the validity of the index: is it smaller than the maximum
	 * value authorized.
	 */
	bool valid();
	/**
	 * Returns the maximum value for the flattened index.
	 */
	long size();

	const Vector<dim,long>& get() const {
		return idx_vector;
	}

private:
	Vector<dim,long> idx_vector;
	Vector<dim,long> n_points_per_dim;
	Vector<dim,long> flatten_weights;
	Vector<dim,long> shift;

	long n_points;

	bool last_idx_fastest;
};

template <int dim>
MultiDimIndex<dim>::MultiDimIndex(
		const Vector<dim,long> &n_points_per_dim,
		const Vector<dim,long> &shift,
		bool last_idx_fastest) :
	n_points_per_dim(n_points_per_dim),
	idx_vector(shift),
	shift(shift),
	last_idx_fastest(last_idx_fastest) {

	if(last_idx_fastest) {
		flatten_weights(dim-1) = 1;
		for(int d=dim-2; d>=0; --d) {
			flatten_weights(d) = flatten_weights(d+1)*n_points_per_dim(d+1);
		}
	}
	else {
		flatten_weights(0) = 1;
		for(int d=1; d<dim; ++d) {
			flatten_weights(d) = flatten_weights(d-1)*n_points_per_dim(d-1);
		}
	}

	n_points=1;
	for(int d=0; d<dim; ++d) {
		n_points *= n_points_per_dim(d);
	}
}

template <int dim>
long& MultiDimIndex<dim>::operator()(long i) {
	return idx_vector(i);
}

template <int dim>
const long& MultiDimIndex<dim>::operator()(long i) const {
	return idx_vector(i);
}

template <int dim>
long MultiDimIndex<dim>::operator()() const {
	return (idx_vector,flatten_weights);
}

template <int dim>
long MultiDimIndex<dim>::operator()(MultiDimIndex<dim> &multi_idx) const {
	return (idx_vector+multi_idx.get(),flatten_weights);
}

template <int dim>
MultiDimIndex<dim>& MultiDimIndex<dim>::operator=(long idx) {

	Assert(
		idx-(flatten_weights,shift)<n_points,
		"Multidim index out of range");
	long rem = idx-(flatten_weights,shift);
	if(last_idx_fastest) {
		for(int d=0; d<dim; d++) {
			idx_vector(d) = shift(d) + floor(rem*1./flatten_weights(d));
			rem = rem % flatten_weights(d);
		}
	}
	else {
		for(int d=dim-1; d>=0; d--) {
			idx_vector(d) = shift(d) + floor(rem*1./flatten_weights(d));
			rem = rem % flatten_weights(d);
		}
	}

	return *this;
}

template <int dim>
MultiDimIndex<dim>& MultiDimIndex<dim>::operator=(
		const MultiDimIndex<dim> &multi_idx) {
	 
	n_points_per_dim = multi_idx.n_points_per_dim;
	flatten_weights = multi_idx.flatten_weights;
	idx_vector = multi_idx.idx_vector;
	n_points = multi_idx.n_points;
	shift = multi_idx.shift;
	last_idx_fastest = multi_idx.last_idx_fastest;

	return *this;
}

template <int dim>
MultiDimIndex<dim>& MultiDimIndex<dim>::operator++() {
	
	int d;
	if(last_idx_fastest) {
		for(d=dim-1; d>0; d--) {
			if(idx_vector(d)-shift(d)<n_points_per_dim(d)-1) {
				idx_vector(d)++;
				break;
			}
			else
				idx_vector(d) = shift(d);
		}
		if(d==0)
			idx_vector(d)++;
	}
	else {
		for(d=0; d<dim-1; d++) {
			if(idx_vector(d)-shift(d)<n_points_per_dim(d)-1) {
				idx_vector(d)++;
				break;
			}
			else
				idx_vector(d) = shift(d);
		}
		if(d==dim-1)
			idx_vector(d)++;
	}

	return *this;
}

template <int dim>
MultiDimIndex<dim> MultiDimIndex<dim>::operator++(int) {
	
	MultiDimIndex<dim> res(*this);
	++(*this);
	return res;
}

template <int dim>
long MultiDimIndex<dim>::size() {
	return n_points;
}

template <int dim>
bool MultiDimIndex<dim>::valid() {
	bool is_valid = true;
	for(int d=0; d<dim; d++) {
		if(idx_vector(d)-shift(d)>=n_points_per_dim(d)) {
			is_valid = false;
			break;
		}
	}
	return is_valid;
}

/**
 * A convenient alias declaration to find the common type of two vectors
 * of different types
 */
template <int dim, typename T1, typename T2> 
using Vec = Vector<dim,typename std::common_type<T1,T2>::type >;

/**
 * Pointwise addition of 2 vectors.
 */
template <int dim, typename T1, typename T2> 
inline Vec<dim,T1,T2> operator+(
		const Vector<dim,T1> &v1, const Vector<dim,T2> &v2) {

	Vec<dim,T1,T2> res;
	for(unsigned int i=0;i<dim; ++i)
		res(i) = v1(i) + v2(i);
	return res;
}

/**
 * Pointwise substraction of 2 vectors.
 */
template <int dim, typename T1, typename T2> 
inline Vec<dim,T1,T2> operator-(
		const Vector<dim,T1> v1, const Vector<dim,T2> &v2) {

	Vec<dim,T1,T2> res;
	for(unsigned int i=0;i<dim; ++i)
		res(i) = v1(i) - v2(i);
	return res;
}

/**
 * Opposite of a vector.
 */
template <int dim, typename T> inline
Vector<dim,T> operator-(const Vector<dim,T> v) {

	Vector<dim,T> res;
	for(unsigned int i=0;i<dim; ++i)
		res(i) = - v(i);
	return res;
}

/**
 * Pointwise multiplication of 2 vectors.
 */
template <int dim, typename T1, typename T2> 
inline Vec<dim,T1,T2> operator*(
		const Vector<dim,T1> v1, const Vector<dim,T2> &v2) {

	Vec<dim,T1,T2> res;
	for(unsigned int i=0;i<dim; ++i)
		res(i) = v1(i) * v2(i);
	return res;
}

/**
 * Pointwise multiplication of a vector and a scalar.
 */
template <int dim, typename T1, typename T2> 
inline Vec<dim,T1,T2> operator*(
		const Vector<dim,T1> v, T2 c) {

	Vec<dim,T1,T2> res;
	for(unsigned int i=0;i<dim; ++i)
		res(i) = v(i) * c;
	return res;
}

/**
 * Pointwise multiplication of a scalar and a vector.
 */
template <int dim, typename T1, typename T2> 
inline Vec<dim,T1,T2> operator*(
		T1 c, const Vector<dim,T2> v) {

	Vec<dim,T1,T2> res;
	for(unsigned int i=0;i<dim; ++i)
		res(i) = v(i) * c;
	return res;
}

/**
 * Pointwise division of 2 vectors.
 */
template <int dim, typename T1, typename T2> 
inline Vec<dim,T1,T2> operator/(
		const Vector<dim,T1> v1, const Vector<dim,T2> &v2) {

	Vec<dim,T1,T2> res;
	for(unsigned int i=0;i<dim; ++i)
		res(i) = v1(i) / v2(i);
	return res;
}

/**
 * Pointwise division of a vector by a scalar.
 */
template <int dim, typename T1, typename T2> 
inline Vec<dim,T1,T2> operator/(
		const Vector<dim,T1> v, T2 c) {

	Vec<dim,T1,T2> res;
	for(unsigned int i=0;i<dim; ++i)
		res(i) = v(i) / c;
	return res;
}

/**
 * Scalar product of 2 vectors.
 */
template <int dim, typename T1, typename T2> 
inline typename std::common_type<T1,T2>::type operator,(
		const Vector<dim,T1> v1, const Vector<dim,T2> &v2) {

	typename std::common_type<T1,T2>::type res = 0;
	for(unsigned int i=0;i<dim; ++i)
		res += v1(i) * v2(i);
	return res;
}

/**
 * Cross-product of 2 three-dimensional vectors.
 */
template <typename T1, typename T2> 
inline Vec<3,T1,T2> operator^(
		const Vector<3,T1> &v1, const Vector<3,T2> &v2) {

	Vec<3,T1,T2> res;
	res(0) = v1(1)*v2(2)-v1(2)*v2(1);
	res(1) = v1(2)*v2(0)-v1(0)*v2(2);
	res(2) = v1(0)*v2(1)-v1(1)*v2(0);
	return res;
}

/**
 * Exponentiation of a Vector by a multidim index.
 */
template <int dim, typename T>
inline T pow(const Vector<dim,T> &v, const MultiDimIndex<dim> &multi_idx) {

	T res = 1.;
	for(int d=0; d<dim; d++) {
		res *= std::pow(v(d), multi_idx(d));
	}
	return res;
}

/**
 * Pretty-printing operator for a vector
 */
template <int dim, typename T> inline
std::ostream& operator<<(std::ostream& os, const Vector<dim,T> &v) {

	os << std::scientific;
	os.precision(8);

	os << std::endl << "┌  ";
	os.fill(' ');
	os.width(18*dim);
	os << " ";
	os << "┐" << std::endl;

	os << "│  ";
	for(unsigned int i=0; i<dim; ++i) {
		os.fill(' ');
		os.width(18);
		os << std::left << v(i);
	}
	os << "│" << std::endl;

	os << "└  ";
	os.fill(' ');
	os.width(18*dim);
	os << " ";
	os << "┘" << std::endl;
	
	return os;
}

/**
 * Parse a Vector from a json object
 */
template <int dim, typename T>
Vector<dim,T> parse_Vector(json j, std::string attribute) {

	std::vector<T> v = j.at(attribute);
	if(v.size()!=dim)
		throw(std::string("Wrong array size for the value of \""+attribute+"\""));

	Vector<dim,int> res = v;
	return res;
}

/**
 * Extract the real part of a complex vector
 */
template <int dim>
Vector<dim,double> real(const Vector<dim,std::complex<double> > &v) {

	Vector<dim,double> w;
	for(int i=0; i<dim; i++)
		w(i) = std::real(v(i));
	return w;
}

/**
 * Return the conjugate of a complex vector
 */
template <int dim>
Vector<dim,std::complex<double> > conj(
		const Vector<dim,std::complex<double> > &v) {

	Vector<dim,std::complex<double> > w;
	for(int i=0; i<dim; i++)
		w(i) = std::conj(v(i));
	return w;
}

/**
 * A convenient alias declaration to find the common type of two
 * matrices of different types
 */
template <int dim1, int dim2, typename T1, typename T2> 
using Mat = Matrix<dim1,dim2,typename std::common_type<T1,T2>::type >;

/**
 * Pointwise addition of two matrices.
 */
template <int dim1, int dim2, typename T1, typename T2>
inline Mat<dim1,dim2,T1,T2> operator+(
		const Matrix<dim1,dim2,T1> &m1, const Matrix<dim1,dim2,T2> &m2) {

	Mat<dim1,dim2,T1,T2> res;
	for(unsigned int i=0; i<dim1; ++i)
		for(unsigned int j=0; j<dim1; ++j)
			res(i,j) = m1(i,j) + m2(i,j);
	return res;
}

/**
 * Pointwise substraction of two matrices.
 */
template <int dim1, int dim2, typename T1, typename T2>
inline Mat<dim1,dim2,T1,T2> operator-(
		const Matrix<dim1,dim2,T1> &m1, const Matrix<dim1,dim2,T2> &m2) {

	Mat<dim1,dim2,T1,T2> res;
	for(unsigned int i=0; i<dim1; ++i)
		for(unsigned int j=0; j<dim1; ++j)
			res(i,j) = m1(i,j) - m2(i,j);
	return res;
}

/**
 * Opposite of a matrix.
 */
template <int dim1, int dim2, typename T>
inline Matrix<dim1,dim2,T> operator-(const Matrix<dim1,dim2,T> m) {
	Matrix<dim1,dim2,T> res;
	for(unsigned int i=0; i<dim1; ++i)
		for(unsigned int j=0; j<dim2; ++j)
			res(i,j) = - m(i,j);
	return res;
}

/**
 * Right matrix-vector product.
 */
template <int dim1, int dim2, typename T1, typename T2>
inline Vec<dim1,T1,T2> operator*(
		const Matrix<dim1,dim2,T1> &m, const Vector<dim2,T2> &v) {

	Vec<dim1,T1,T2> res;
	for(unsigned int i=0; i<dim1; ++i)
		for(unsigned int j=0; j<dim2; ++j)
			res(i) += m(i,j) * v(j);
	return res;
}

/**
 * Left matrix-vector product.
 */
template <int dim1, int dim2, typename T1, typename T2>
inline Vec<dim2,T1,T2> operator*(
		const Vector<dim1,T1> &v, const Matrix<dim1,dim2,T2> &m) {

	Vec<dim2,T1,T2> res;
	for(unsigned int i=0; i<dim2; ++i)
		for(unsigned int j=0; j<dim1; ++j)
			res(i) += v(j) * m(j,i);
	return res;
}

/**
 * Pointwise multiplication of a matrix and a scalar.
 */
template <int dim1, int dim2, typename T1, typename T2>
inline Mat<dim1,dim2,T1,T2> operator*(
		const Matrix<dim1,dim2,T1> &m, T2 c) {

	Mat<dim1,dim2,T1,T2> res;
	for(unsigned int i=0; i<dim1; ++i)
		for(unsigned int j=0; j<dim2; ++j)
			res(i,j) = m(i,j) * c;
	return res;
}

/**
 * Pointwise multiplication of a scalar and a matrix.
 */
template <int dim1, int dim2, typename T1, typename T2>
inline Mat<dim1,dim2,T1,T2> operator*(
		T1 c, const Matrix<dim1,dim2,T2> &m) {

	Mat<dim1,dim2,T1,T2> res;
	for(unsigned int i=0; i<dim1; ++i)
		for(unsigned int j=0; j<dim2; ++j)
			res(i,j) = m(i,j) * c;
	return res;
}

/**
 * Pretty-printing operator for a matrix.
 */
template <int dim1, int dim2, typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<dim1,dim2,T> &m) {

	os << std::scientific;
	os.precision(2);

	os << std::endl << "┌  ";
	os.fill(' ');
	os.width(12*dim2);
	os << " ";
	os << "┐" << std::endl;

	for(unsigned int i=0; i<dim1; ++i) {
		os << "│  ";
		for(unsigned int j=0; j<dim2; ++j) {
			os.fill(' ');
			os.width(12);
			os << std::left << m(i,j);
		}
		os << "│" << std::endl;
	}

	os << "└  ";
	os.fill(' ');
	os.width(12*dim2);
	os << " ";
	os << "┘" << std::endl;
	
	return os;
}

/**
 * @}
 */

///////////////////////////////////////////////
// Inline implementation of member operators //
///////////////////////////////////////////////

template <int dim, typename T1>
inline Vector<dim,T1>::Vector() {

	for(unsigned int i=0;i<dim; ++i)
		data[i] = 0;
}

template <int dim, typename T1> 
template <typename T2>
inline Vector<dim,T1>::Vector(T2 c) {

	for(unsigned int i=0;i<dim; ++i)
		data[i] = c;
}

template <int dim, typename T1>
template <typename T2>
inline Vector<dim,T1>::Vector(const Vector<dim,T2> &v) {

	for(unsigned int i=0; i<dim; ++i)
		data[i] = v(i);
}

template <int dim, typename T1> 
template <typename T2>
inline Vector<dim,T1>::Vector(const T2 (&v)[dim]) {

	for(unsigned int i=0; i<dim; ++i)
		data[i] = v[i];
}

template <int dim, typename T1> 
template <typename T2>
inline Vector<dim,T1>::Vector(const std::vector<T2> &v) {

	if(dim!=v.size())
		throw(std::string(
			"You are trying to initialize a Vector with an std::vector of "
			"the wrong size"));
	for(unsigned int i=0;i<dim; ++i)
		data[i] = v[i];
}

template <int dim, typename T1>
inline double Vector<dim,T1>::norm() const {
	
	double res = 0;
	for(unsigned int i=0; i<dim; ++i)
		res += std::abs(data[i])*std::abs(data[i]);
	return std::sqrt(res);
}

template <int dim, typename T1> 
template <typename T2>
inline Vector<dim,T1>& Vector<dim,T1>::operator=(const Vector<dim,T2> &v) {

	for(unsigned int i=0;i<dim; ++i)
		data[i] = v(i);
	return *(this);
}

template <int dim, typename T1>
template <typename T2>
inline Vector<dim,T1>& Vector<dim,T1>::operator=(const T2 (&v)[dim]) {

	for(unsigned int i=0;i<dim; ++i)
		data[i] = v[i];
	return *(this);
}

template <int dim, typename T1>
template <typename T2>
inline Vector<dim,T1>& Vector<dim,T1>::operator=(const std::vector<T2> &v) {

	if(dim!=v.size())
		throw(std::string(
			"You are trying to initialize a Vector with an std::vector of "
			"the wrong size"));
	for(unsigned int i=0;i<dim; ++i)
		data[i] = v[i];
	return *(this);
}

template <int dim, typename T1>
inline Vector<dim,T1>& Vector<dim,T1>::operator=(T1 c) {

	for(unsigned int i=0;i<dim; ++i)
		data[i] = c;
	return *(this);
}

template <int dim, typename T1>
template <typename T2>
inline Vector<dim,T1>& Vector<dim,T1>::operator+=(const Vector<dim,T2> &v) {

	for(unsigned int i=0;i<dim; ++i)
		data[i] += v(i);
	return *(this);
}

template <int dim, typename T1>
template <typename T2>
inline Vector<dim,T1>& Vector<dim,T1>::operator-=(const Vector<dim,T2> &v) {

	for(unsigned int i=0;i<dim; ++i)
		data[i] -= v(i);
	return *(this);
}

template <int dim, typename T1>
template <typename T2>
inline Vector<dim,T1>& Vector<dim,T1>::operator*=(T2 c) {

	for(unsigned int i=0;i<dim; ++i)
		data[i] *= c;
	return *(this);
}

template <int dim, typename T1>
template <typename T2>
inline Vector<dim,T1>& Vector<dim,T1>::operator/=(T2 c) {

	for(unsigned int i=0;i<dim; ++i)
		data[i] /= c;
	return *(this);
}

template <int dim, typename T1>
inline T1& Vector<dim,T1>::operator()(unsigned int i) {

	Assert(i<dim,
		"The desired index is greater than the dimension of the vector");
	return data[i];
}

template <int dim, typename T1>
inline const T1& Vector<dim,T1>::operator()(unsigned int i) const {

	Assert(i<dim,
		"The desired index is greater than the dimension of the vector");
	return data[i];
}

template <int dim1, int dim2, typename T1>
inline Matrix<dim1,dim2,T1>::Matrix() {

	for(unsigned int i=0;i<dim1; ++i)
		for(unsigned int j=0;j<dim2; ++j)
			data[i+dim1*j] = 0;
}

template <int dim1, int dim2, typename T1>
inline Matrix<dim1,dim2,T1>::Matrix(T1 c) {

	for(unsigned int i=0;i<dim1; ++i)
		for(unsigned int j=0;j<dim2; ++j)
			data[i+dim1*j] = c;
}

template <int dim1, int dim2, typename T1>
template <typename T2>
inline Matrix<dim1,dim2,T1>::Matrix(const Matrix<dim1,dim2,T2> &m) {

	for(unsigned int i=0;i<dim1; ++i)
		for(unsigned int j=0;j<dim2; ++j)
			data[i+dim1*j] = m(i,j);
}

template <int dim1, int dim2, typename T1>
template <typename T2>
inline Matrix<dim1,dim2,T1>::Matrix(const T2 (&m)[dim1*dim2]) {

	for(unsigned int i=0;i<dim1*dim2; ++i)
		data[i] = m[i];
}

template <int dim1, int dim2, typename T1>
inline double Matrix<dim1,dim2,T1>::norm() const {
	
	double res = 0;
	for(unsigned int i=0; i<dim1*dim2; ++i)
		res += std::abs(data[i])*std::abs(data[i]);
	return std::sqrt(res);
}

template <int dim1, int dim2, typename T1>
template <typename T2>
inline Matrix<dim1,dim2,T1>& Matrix<dim1,dim2,T1>::operator=(
		const Matrix<dim1,dim2,T2> &m) {

	for(unsigned int i=0;i<dim1; ++i)
		for(unsigned int j=0;j<dim2; ++j)
			data[i+dim1*j] = m(i,j);
	return *(this);
}

template <int dim1, int dim2, typename T1>
inline Matrix<dim1,dim2,T1>& Matrix<dim1,dim2,T1>::operator=(T1 c) {

	for(unsigned int i=0;i<dim1; ++i)
		for(unsigned int j=0;j<dim2; ++j)
			data[i+dim1*j] = c;
	return *(this);
}

template <int dim1, int dim2, typename T1>
template <typename T2>
inline Matrix<dim1,dim2,T1>& Matrix<dim1,dim2,T1>::operator+=(
		const Matrix<dim1,dim2,T2> &m) {

	for(unsigned int i=0;i<dim1; ++i)
		for(unsigned int j=0;j<dim2; ++j)
			data[i+dim1*j] += m(i,j);
	return *(this);
}

template <int dim1, int dim2, typename T1>
template <typename T2>
inline Matrix<dim1,dim2,T1>& Matrix<dim1,dim2,T1>::operator-=(
		const Matrix<dim1,dim2,T2> &m) {

	for(unsigned int i=0;i<dim1; ++i)
		for(unsigned int j=0;j<dim2; ++j)
			data[i+dim1*j] -= m(i,j);
	return *(this);
}

template <int dim1, int dim2, typename T1>
template <typename T2>
inline Matrix<dim1,dim2,T1>& Matrix<dim1,dim2,T1>::operator*=(T2 c) {
	for(unsigned int i=0;i<dim1; ++i)
		for(unsigned int j=0;j<dim2; ++j)
			data[i+dim1*j] *= c;
	return *(this);
}

template <int dim1, int dim2, typename T1>
inline T1& Matrix<dim1,dim2,T1>::operator()(unsigned int i, unsigned int j) {

	Assert(i<dim1 && j<dim2,
		"Out of bound matrix indices request.");
	return data[i+dim1*j];
}

template <int dim1, int dim2, typename T1>
inline const T1& Matrix<dim1,dim2,T1>::operator()(
		unsigned int i, unsigned int j) const {

	Assert(i<dim1 && j<dim2,
		"Out of bound matrix indices request.");
	return data[i+dim1*j];
}

namespace Basis {
	const Vector<3,double> ex({1,0,0});
	const Vector<3,double> ey({0,1,0});
	const Vector<3,double> ez({0,0,1});
};

template <int dim>
class Identity : public Matrix<dim,dim,double> {
public:
	Identity() : Matrix<dim,dim,double>::Matrix(0) {
		for(int i=0; i<dim; i++)
			this->data[i*(dim+1)] = 1.;
	}
};

template <int dim>
class BasisVector : public Vector<dim,double> {
public:
	BasisVector(int i) : Vector<dim,double>(0.) {
		if(i<0 || i>dim-1)
			throw(std::string("Error: wrong index for the basis vector."));
		this->data[i] = 1;
	}
};

///////////////////////////////////////////////////
// Ugly trick to have a less shitty std::complex //
///////////////////////////////////////////////////
template <typename T>
struct identity_t { typedef T type; };

#define COMPLEX_OPS(OP)                                                                      \
	template <typename _Tp>                                                                  \
	std::complex<_Tp>                                                                        \
	operator OP(std::complex<_Tp> lhs, const typename identity_t<_Tp>::type & rhs) {         \
		return lhs OP rhs;                                                                   \
	}                                                                                        \
	template <typename _Tp>                                                                  \
	std::complex<_Tp>                                                                        \
	operator OP(const typename identity_t<_Tp>::type & lhs, const std::complex<_Tp> & rhs) { \
		return lhs OP rhs;                                                                   \
	}

COMPLEX_OPS(+)
COMPLEX_OPS(-)
COMPLEX_OPS(*)
COMPLEX_OPS(/)
#undef COMPLEX_OPS

#endif
