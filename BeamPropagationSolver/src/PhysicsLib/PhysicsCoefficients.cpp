#include <muParser.h>

#include "PhysicsCoefficients.h"

PhysicsCoefficients::PhysicsCoefficients(
		const RootSettings &settings,
		const CartesianMesh &mesh) :
	coefs_settings(settings.physics.coefs) {}

double PhysicsCoefficients::get_n1(unsigned int domain_id,double wavelength) const {

	try {
		mu::Parser p;
		p.DefineConst("lambda", wavelength);
		p.SetExpr(coefs_settings.nsample_expressions().at(3*domain_id));
		return p.Eval();
	}
	catch(mu::Parser::exception_type &e) {
		throw std::string("Wrong expression for n1");
	}
}

double PhysicsCoefficients::get_n2(unsigned int domain_id,double wavelength) const {

	try {
		mu::Parser p;
		p.DefineConst("lambda", wavelength);
		p.SetExpr(coefs_settings.nsample_expressions().at(3*domain_id+1));
		return p.Eval();
	}
	catch(mu::Parser::exception_type &e) {
		throw std::string("Wrong expression for n2");
	}
}

double PhysicsCoefficients::get_n3(unsigned int domain_id,double wavelength) const {

	try {
		mu::Parser p;
		p.DefineConst("lambda", wavelength);
		p.SetExpr(coefs_settings.nsample_expressions().at(3*domain_id+2));
		return p.Eval();
	}
	catch(mu::Parser::exception_type &e) {
		throw std::string("Wrong expression for n3");
	}
}

double PhysicsCoefficients::get_niso_up(unsigned int layer_id, double wavelength) const {

	try {
		mu::Parser p;
		p.DefineConst("lambda", wavelength);
		p.SetExpr(coefs_settings.niso_up_expressions().at(layer_id));
		return p.Eval();
	}
	catch(mu::Parser::exception_type &e) {
		throw std::string("Wrong expression for niso_up");
	}
}

double PhysicsCoefficients::get_niso_lo(unsigned int layer_id, double wavelength) const {

	try {
		mu::Parser p;
		p.DefineConst("lambda", wavelength);
		p.SetExpr(coefs_settings.niso_lo_expressions().at(layer_id));
		return p.Eval();
	}
	catch(mu::Parser::exception_type &e) {
		throw std::string("Wrong expression for niso_lo");
	}
}

double PhysicsCoefficients::get_nin(double wavelength) const {

	try {
		mu::Parser p;
		p.DefineConst("lambda", wavelength);
		p.SetExpr(coefs_settings.nin_expression);
		return p.Eval();
	}
	catch(mu::Parser::exception_type &e) {
		throw std::string("Wrong expression for nin");
	}
}

double PhysicsCoefficients::get_nout(double wavelength) const {

	try {
		mu::Parser p;
		p.DefineConst("lambda", wavelength);
		p.SetExpr(coefs_settings.nout_expression);
		return p.Eval();
	}
	catch(mu::Parser::exception_type &e) {
		throw std::string("Wrong expression for nout");
	}
}
