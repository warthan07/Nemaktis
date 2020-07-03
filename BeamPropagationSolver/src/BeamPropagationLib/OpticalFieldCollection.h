#ifndef OPTICALFIELDCOLLECTION_H
#define OPTICALFIELDCOLLECTION_H

#include <complex>
#include <utility>

#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkPointData.h>

#include "CartesianMesh.h"
#include "PhysicsCoefficients.h"
#include "OpticalField.h"

class ScreenOpticalFieldCollection {
public:
	ScreenOpticalFieldCollection(
		const CartesianMesh &mesh,
		const PhysicsCoefficients &coefs);

	TransverseOpticalField& operator()(int wave_idx, int q_idx, int pol_idx) {
		return fields[pol_idx+2*(q_idx+q_vals.size()*wave_idx)];
	}

	const CartesianMesh mesh;

private:
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
	 * polarisation.
	 */
	std::vector<TransverseOpticalField> fields;
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
