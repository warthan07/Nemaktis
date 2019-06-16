#include <muParser.h>

#include "PhysicsCoefficients.h"

PhysicsCoefficients::PhysicsCoefficients(
		const RootSettings &settings,
		const CartesianMesh &mesh) :
	physics_settings(settings.physics),
    _wavelengths(settings.physics.coefs.wavelengths()) {

	ne_expression = physics_settings.coefs.ne_expression;
	no_expression = physics_settings.coefs.no_expression;
	nhost_expression = physics_settings.coefs.nhost_expression;

	_mesh_volume = 
		mesh.delta_x*(mesh.Nx-1)*
		mesh.delta_y*(mesh.Ny-1)*
		mesh.delta_z*(mesh.Nz-1);
	_mesh_thickness = mesh.delta_x*(mesh.Nx-1);
	_z_origin = -mesh.delta_z*(mesh.Nz-1)/2.;
}

double PhysicsCoefficients::get_ne(double wavelength) const {

	try {
		mu::Parser p;
		p.DefineConst("lambda", wavelength);
		p.SetExpr(ne_expression);
		return p.Eval();
	}
	catch(mu::Parser::exception_type &e) {
		throw std::string("Wrong expression for ne");
	}
}

double PhysicsCoefficients::get_no(double wavelength) const {

	try {
		mu::Parser p;
		p.DefineConst("lambda", wavelength);
		p.SetExpr(no_expression);
		return p.Eval();
	}
	catch(mu::Parser::exception_type &e) {
		throw std::string("Wrong expression for no");
	}
}

double PhysicsCoefficients::get_nhost(double wavelength) const {

	try {
		mu::Parser p;
		p.DefineConst("lambda", wavelength);
		p.SetExpr(nhost_expression);
		return p.Eval();
	}
	catch(mu::Parser::exception_type &e) {
		throw std::string("Wrong expression for nhost");
	}
}
