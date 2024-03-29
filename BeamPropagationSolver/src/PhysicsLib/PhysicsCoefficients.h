/** 
 * @file PhysicsCoefficients.h
 */

#ifndef LCCOEFFICIENTS_H
#define LCCOEFFICIENTS_H

#include <memory>
#include <map>

#include "RootSettings.h"
#include "CartesianMesh.h"

#define PI 3.1415926535897932

class PhysicsCoefficients {
public:
	PhysicsCoefficients(
		const RootSettings &settings,
		const CartesianMesh &mesh);

	double get_ne(double wavelength) const;
	double get_ne_imag(double wavelength) const;
	double get_no(double wavelength) const;
	double get_nhost(double wavelength) const;
	double get_nin(double wavelength) const;
	double get_nout(double wavelength) const;

	double mesh_volume() const {
		return _mesh_volume;
	}
	double mesh_thickness() const {
		return _mesh_thickness;
	}
	double z_origin() const {
		return _z_origin;
	}
	std::vector<double> wavelengths() const {
		return _wavelengths;
	}
	std::vector<std::pair<double,double> > q_vals() const {
		return _q_vals;
	}

private:
	std::string ne_expression;
	std::string ne_imag_expression;
	std::string no_expression;
	std::string nhost_expression;
	std::string nin_expression;
	std::string nout_expression;

	const PhysicsSettings &physics_settings;

	/**
	 * Volume of the LC layer.
	 */
	double _mesh_volume;
	/**
	 * Thickness of the LC layer.
	 */
	double _mesh_thickness;
	/**
	 * Z-origin for the mesh
	 */
	double _z_origin;

	/**
	 * Array containing all the wavelengths in the light spectrum.
	 */
	std::vector<double> _wavelengths;
	/**
	 * Array containing x and y components of the incoming wavectors.
	 */
	std::vector<std::pair<double,double> > _q_vals;
};

#endif
