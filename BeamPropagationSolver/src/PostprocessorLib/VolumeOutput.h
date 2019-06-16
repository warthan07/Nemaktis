#ifndef VOLUMEOUTPUT_H
#define VOLUMEOUTPUT_H

#include "PostprocessorBase.h"

class VolumeOutput : public PostprocessorBase {
public:
	VolumeOutput(
		const RootSettings &settings,
		const PhysicsCoefficients &coefs);

	virtual void apply(
		VectorField<double> &lc_sol,
		std::vector<VectorField<std::complex<double> > > (&bpm_sol)[2]);

private:
	/**
	 * Array containing all the wavelengths in the light spectrum.
	 */
	std::vector<double> wavelengths;
};
#endif
