#ifndef ALGORITHMSETTINGS_H
#define ALGORITHMSETTINGS_H

#include "BaseSettings.h"

enum class BoundaryType {
	Transparent, Periodic
};

class GeneralSettings : public BaseSettings {
public:
	GeneralSettings(const nlohmann::json &j);

	/**
	 * Type of optical axis field for initializing the permittivity tensor field. Should be
	 * None (only isotropic media), Director (isotropic or uniaxial media), QTensor
	 * (isotropic or uniaxial media, with possible biaxiality near defects' cores) or
	 * DualDirector (isotropic, uniaxial or biaxial media)
	 */
	const std::string optical_axis_field_type;
	/**
	 * Path name for the results folder, which will be created if inexistent.
	 */
	const std::string results_folder_name;
};

class BPMSettings : public BaseSettings {
public:
	BPMSettings(const nlohmann::json &j);

	/**
	 * Number of woodbury iterations for the iterative XY corrections in the diffractions
	 * operators. For most systems, a small value such as 2 is enough.
	 */
	const int N_woodbury_steps;
	/**
	 * Number of sub-integration steps for the evolution operators. A number higher than 1
	 * allows a better evaluation of the diffraction operator, but will obviously take more
	 * time.
	 */
	const int Nz_substeps;
};

class AlgorithmSettings : public BaseSettings {
public:
	AlgorithmSettings(const nlohmann::json &j);

	/**
	 * General settings for the software
	 */
	GeneralSettings general;
	/**
	 * Settings for the beam propagation algorithm
	 */
	BPMSettings bpm;
};

#endif
