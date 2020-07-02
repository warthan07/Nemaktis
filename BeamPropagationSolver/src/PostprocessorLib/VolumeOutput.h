#ifndef VOLUMEOUTPUT_H
#define VOLUMEOUTPUT_H

#include "PostprocessorBase.h"
#include "OpticalFieldCollection.h"

class VolumeOutput : public PostprocessorBase {
public:
	VolumeOutput(
		const RootSettings &settings,
		const PhysicsCoefficients &coefs);

	void apply(BulkOpticalFieldCollection &bulk_optical_fields);
};
#endif
