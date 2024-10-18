#ifndef OPTICALFIELDCOLLECTION_H
#define OPTICALFIELDCOLLECTION_H

#include <complex>
#include <memory>
#include <utility>

#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkPointData.h>

#include <fftw3.h>

#include "CartesianMesh.h"
#include "PhysicsCoefficients.h"
#include "OpticalField.h"

class ScreenOpticalFieldCollection {
public:
	ScreenOpticalFieldCollection(
		const CartesianMesh &mesh,
		const PhysicsCoefficients &coefs);

	ScreenOpticalFieldCollection(
		const CartesianMesh &mesh,
		const PhysicsCoefficients &coefs,
		std::complex<double>* user_vals,
		unsigned int n_user_vals);

	~ScreenOpticalFieldCollection() {
		if(own_field_vals)
			fftw_free(field_vals);
		fftw_free(fft_field_vals);
	}

	TransverseOpticalField& operator()(int wave_idx, int q_idx, int pol_idx) {
		return *fields[pol_idx+2*(q_idx+q_vals.size()*wave_idx)];
	}
	TransverseOpticalField& fft(int wave_idx, int q_idx, int pol_idx) {
		return *fft_fields[pol_idx+2*(q_idx+q_vals.size()*wave_idx)];
	}

	void apply_forward_fft_plan(int wave_idx) {
		fftw_execute_dft(
			forward_plan, &field_vals[wave_idx*stride_plans], &fft_field_vals[wave_idx*stride_plans]);
	}
	void apply_backward_fft_plan(int wave_idx) {
		fftw_execute_dft(
			backward_plan, &fft_field_vals[wave_idx*stride_plans], &field_vals[wave_idx*stride_plans]);
	}

	const CartesianMesh mesh;

private:
	/**
	 * Initialize all transverse optical fields and fft plans
	 */
	void init_fields_and_plans();

	/**
	 * Array containing all the wavelengths in the light spectrum.
	 */
	std::vector<double> wavelengths;
	/**
	 * Array containing the incoming wavectors.
	 */
	std::vector<std::pair<double,double> > q_vals;

	/**
	 * Vectors of transverse optical fields for each wavelength, wavevector, and
	 * polarisation. The underlying data vectors are aligned and contiguous in memory and
	 * pointed by the private attribute field_vals and fft_field_vals that can be used with fftw.
	 */
	std::vector<std::shared_ptr<TransverseOpticalField> > fields, fft_fields;

	bool own_field_vals;
	int stride_fields, stride_plans;
	fftw_complex *field_vals, *fft_field_vals;
	fftw_plan forward_plan, backward_plan;
};

class BulkOpticalFieldCollection {
public:
	BulkOpticalFieldCollection(
		const CartesianMesh &mesh,
		const PhysicsCoefficients &coefs);

	void set_field_val(
		int wave_idx, int q_idx, int pol_idx, int comp, int mesh_idx,
		std::complex<double> val);

	vtkSmartPointer<vtkImageData> get_vti_data() const {
		return vti_data;
	}
	const std::vector<double> get_wavelengths() const {
		return wavelengths;
	}
	const std::vector<std::pair<double,double> > get_q_vals() const {
		return q_vals;
	}

private:
	/**
	 * Mesh dimensions
	 */
	int Nx, Ny, Nz;

	/**
	 * VTK structure storing all the arrays of optical fields.
	 */
	vtkSmartPointer<vtkImageData> vti_data;
	/**
	 * Vectors of VTK array for the real part of the optical fields.
	 */
	std::vector<vtkSmartPointer<vtkDoubleArray> > real_data_arrays;
	/**
	 * Vectors of VTK array for the imaginary part of the optical fields.
	 */
	std::vector<vtkSmartPointer<vtkDoubleArray> > imag_data_arrays;

	/**
	 * Array containing all the wavelengths in the light spectrum.
	 */
	std::vector<double> wavelengths;
	/**
	 * Array containing the incoming wavectors.
	 */
	std::vector<std::pair<double,double> > q_vals;
};

#endif
