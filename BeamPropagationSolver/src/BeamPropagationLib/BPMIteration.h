#ifndef BPMITERATION_H
#define BPMITERATION_H

#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

#include "BeamProfile.h"
#include "VectorField.h"
#include "OpticalFieldCollection.h"

class BPMIteration {
public:
	BPMIteration(
		const VectorField<double> &lc_sol,
		ScreenOpticalFieldCollection &screen_optical_fields,
		const PhysicsCoefficients &coefs,
		const RootSettings &settings);

	void propagate_fields();

	std::shared_ptr<BulkOpticalFieldCollection> get_bulk_optical_fields() {
		return bulk_optical_fields;
	}

private:
	/**
	 * Type of the beam profile function.
	 */
	BeamProfileType beam_profile_type;

	/**
	 * Reference to the PhysicsCoefficients object of this run.
	 */
	const PhysicsCoefficients &coefs;

	/**
	 * Reference to the RootSettings object of this run.
	 */
	const RootSettings &settings;

	/**
	 * Mesh spacings.
	 */
	double delta_x, delta_y, delta_z;
	/**
	 * Number of points in each space direction for the global mesh
	 */
	int Nx, Ny, Nz;
	/**
	 * Number of evolution sub-steps per z-slab
	 */
	unsigned int Nz_substeps;

	/**
	 * Beam waist (in case of a gaussian input
	 */
	double waist;

	/**
	 * Shortcut reference to the LC orientational field (n or q).
	 */
	const VectorField<double> &lc_field;

	/**
	 * Reference to the object storing the values of the optical fields for the current
	 * transverse plane.
	 */
	ScreenOpticalFieldCollection &screen_optical_fields;

	/**
	 * Do we need to export bulk values of the optical field?
	 */
	bool bulk_output;
	/**
	 * Object storing the bulk values of the optical field, if bulk_output is true.
	 */
	std::shared_ptr<BulkOpticalFieldCollection> bulk_optical_fields;

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
