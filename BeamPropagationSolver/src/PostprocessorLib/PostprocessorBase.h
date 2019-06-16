/**
 * @file PostprocessorBase.h
 */

#ifndef POSTPROCESSORBASE_H
#define POSTPROCESSORBASE_H

#include <complex>

#include <vtkSmartPointer.h>
#include <vtkImageData.h>

#include "PhysicsCoefficients.h"
#include "PostprocessorSettings.h"
#include "VectorField.h"

class PostprocessorBase {
public:
	/**
	 * Constructor.
	 */
	PostprocessorBase(
		const BasePostprocessorSettings &pp_settings,
		const std::string &results_dirname,
		const PhysicsCoefficients &coefs);

	/**
	 * Postprocessing operator, defined here as a pure virtual method.
	 */
	virtual void apply(
		VectorField<double> &lc_sol,
		std::vector<VectorField<std::complex<double> > > (&bpm_sol)[2]) = 0;

protected:
	/**
	 * A method which format the main output path+filename from
	 * base_name and the current values of the variable coefs.
	 */
	std::string format_global_path(std::string suffix = "");

	/**
	 * The base name for the exported file or directory.
	 */
	const std::string base_name;
	/**
	 * The path to where the results will be exported.
	 */
	const std::string results_dirname;

	/**
	 * Reference  to the PhysicsCoefficients object of the problem.
	 */
	const PhysicsCoefficients& coefs;
};

#endif
