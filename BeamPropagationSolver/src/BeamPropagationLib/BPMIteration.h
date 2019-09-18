#ifndef BPMITERATION_H
#define BPMITERATION_H

#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

#include "BeamProfile.h"
#include "VectorField.h"

class BPMIteration {
public:
	BPMIteration(
		const VectorField<double> &lc_sol,
		std::vector<VectorField<std::complex<double> > > (&bpm_sol)[2],
		const PhysicsCoefficients &coefs,
		const RootSettings &settings);

	void update_optical_field();

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
	 * Reference to the BPMSettings object of this run.
	 */
	const BPMSettings &bpm_settings;

	/**
	 * Mesh spacings.
	 */
	double delta_x, delta_y, delta_z;
	/**
	 * Number of points in each space direction for the global mesh
	 * (LC+plates)
	 */
	int Nx, Ny, Nz;

	/**
	 * Beam waist (in case of a gaussian input
	 */
	double waist;

	/**
	 * Shortcut reference to the LC solution.
	 */
	const VectorField<double> &lc_sol;
	/**
	 * Shortcut reference to the BPM solution.
	 */
	std::vector<VectorField<std::complex<double> > > (&bpm_sol)[2];

	/**
	 * Array containing all the wavelengths in the light spectrum.
	 */
	std::vector<double> wavelengths;
};

#endif
