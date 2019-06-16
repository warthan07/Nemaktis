#ifndef ROOTSETTINGS_H
#define ROOTSETTINGS_H

#include "BaseSettings.h"
#include "AlgorithmSettings.h"
#include "PhysicsSettings.h"
#include "PostprocessorSettings.h"

class RootSettings : public BaseSettings {
	public:
		RootSettings(const nlohmann::json &j);

		AlgorithmSettings algorithm;
		PhysicsSettings physics;
		PostprocessorSettings postprocessor;
};

void runFromSettings(RootSettings &settings);

#endif
