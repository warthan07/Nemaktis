#ifndef HCINVERSESOLVER_H
#define HCINVERSESOLVER_H

#include "CubicInterpolatedMapping.h"

template <int dim>
using DataPair = std::pair<
	Vector<dim,double>,
	std::vector<Vector<dim,double> > >;

template <int dim>
class HCInverseSolver {
public:
	HCInverseSolver(CubicInterpolatedMapping<dim,dim,double> map, double typical_length, json settings);

	void find_inverses(
		Vector<dim,double> &y, std::vector<Vector<dim,double> > &sols); 

	void find_inverses(
		DataPair<dim>& output, DataPair<dim>& input);
	
protected:
	bool predictor_step(
		std::vector<Vector<dim,double> > &sols);
	bool corrector_step(
		std::vector<Vector<dim,double> > &sols);
	bool newton_update();

	virtual bool update_residual() = 0;
	virtual bool update_jacobian_and_residual() = 0;
	bool update_reduced_jacobian_and_residual();
	bool compute_root_func(const Vector<dim,double> &pos, Vector<dim,double> &res);

	double newton_tol;
	int opt_iter;
	double max_step;
	double max_arclength;

	double jac_perp_det;
	Matrix<dim,dim+1,double> jac;
	Matrix<dim,dim,double> jac_perp, jac_perp_inv;
	Matrix<dim+1,dim,double> nullspace_mat;

	Vector<3,double> e1,e2;
	Vector<dim+1,double> corrector_update, predictor_update;
	Vector<dim+1,double> hc_pos, prev_hc_pos, old_hc_pos;

	Vector<dim,double> res;

	double arclength_step, boundary_step;
	bool record_sol, loop_test, loop_detected;

	CubicInterpolatedMapping<dim,dim,double> map;
	Vector<dim,double> y;
	Vector<dim,double> x;
	Vector<dim,double> x0;
	double L;

	Vector<dim,double> map_val, map_val_x0;
	Matrix<dim,dim,double> map_grad;
};

template <int dim>
class FPHCInverseSolver : public HCInverseSolver<dim> {
public:
	FPHCInverseSolver(
			CubicInterpolatedMapping<dim,dim,double> map, double typical_length, json settings) :
		HCInverseSolver<dim>::HCInverseSolver(map, typical_length, settings) {};

protected:
	virtual bool update_residual();
	virtual bool update_jacobian_and_residual();
};

template <int dim>
class NHCInverseSolver : public HCInverseSolver<dim> {
public:
	NHCInverseSolver(
			CubicInterpolatedMapping<dim,dim,double> map, double typical_length, json settings) :
		HCInverseSolver<dim>::HCInverseSolver(map, typical_length, settings) {};

protected:
	virtual bool update_residual();
	virtual bool update_jacobian_and_residual();
};

#endif
