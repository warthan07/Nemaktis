#include <fstream>

#include "HCInverseSolver.h"

template <int dim>
HCInverseSolver<dim>::HCInverseSolver(
		CubicInterpolatedMapping<dim,dim,double> map, double typical_length, json j) :
	map(map),
	L(typical_length) {

	newton_tol = j.at("Newton tolerance");
	opt_iter = j.at("Optimum newton step number");
	max_step = L * j.at("Max arclength step").get<double>();
	max_arclength = L * j.at("Max arclength").get<double>();
}

template <int dim>
void HCInverseSolver<dim>::find_inverses(
		Vector<dim,double> &_y, std::vector<Vector<dim,double> > &sols) {

	y = _y;
	x0 = y;

	if(!map.get_value(x0,map_val_x0)) {
		map.get_def_domain()->compute_projection(y, x0);
		if(!map.get_value(x0,map_val_x0))
			throw(std::string(
				"Error: impossible to inverse the screen mapping because one of the \n"
				"screen point is outside the definition domain of the light source. \n"
				"You must specify a bigger light source."));
	}

	sols.clear();

	// We explore the forward branch
	x = x0;
	arclength_step = max_step;
	loop_test = true;
	loop_detected = false;

	for(int i=0; i<dim; i++)
		old_hc_pos(i) = x(i);
	old_hc_pos(dim) = -arclength_step; // force (dtau/ds)>0

	hc_pos = old_hc_pos;
	hc_pos(dim) = 0;

	double arclength = 0;
	while(predictor_step(sols) && corrector_step(sols)) {
		arclength += (hc_pos-old_hc_pos).norm();
		if(sols.size() > 2 || arclength > max_arclength)
			break;
	}

	// We explore the backward branch
	if(!loop_detected) {
		x = x0;
		arclength_step = max_step;
		loop_test = true;
		loop_detected = false;

		for(int i=0; i<dim; i++)
			old_hc_pos(i) = x(i);
		old_hc_pos(dim) = arclength_step; // force (dtau/ds)<0

		hc_pos = old_hc_pos;
		hc_pos(dim) = 0;

		arclength = 0;
		while(predictor_step(sols) && corrector_step(sols)) {
			arclength += (hc_pos-old_hc_pos).norm();
			if(sols.size()>2 || arclength > max_arclength)
				break;
		}
	}
}

template <int dim>
void HCInverseSolver<dim>::find_inverses(
		DataPair<dim>& output, DataPair<dim>& input) {

	y = output.first;

	Matrix<dim,dim,double> grad, grad_inv;
	Vector<dim,double> shift;
	std::size_t prev_size;
	
	output.second.clear();

	for(auto prev_sol : input.second) {
		x0 = prev_sol;
		if(!map.get_value(x0,map_val_x0))
			throw(std::string(
				"Error: impossible to inverse the screen mapping because one of the \n"
				"starting points of the HC algorithm is outside the definition\n"
				"domain of the light source. You should probably specify a larger\n"
				"light source."));

		// We explore the forward branch
		x = x0;
		arclength_step = max_step;
		loop_test = true;
		loop_detected = false;

		for(int i=0; i<dim; i++)
			old_hc_pos(i) = x(i);
		old_hc_pos(dim) = -arclength_step; // force (dtau/ds)>0

		hc_pos = old_hc_pos;
		hc_pos(dim) = 0;

		double arclength = 0;
		prev_size = output.second.size();
		while(predictor_step(output.second) && corrector_step(output.second)) {
			arclength += (hc_pos-old_hc_pos).norm();
			if(output.second.size() > prev_size || arclength > max_arclength)
				break;
		}
	}
}

template <int dim>
bool HCInverseSolver<dim>::predictor_step(
		std::vector<Vector<dim,double> > &sols) {

	if(!update_jacobian_and_residual()) {
		return false;
	}

	// We compute the hc_update in the extended homotopy space (x,tau)
	// and the nullspace matrix containing a basis for the hyperspace
	// orthogonal to the predictor_update
	if(dim==1) {
		predictor_update(0) = -jac(0,1);
		predictor_update(1) = jac(0,0);

		double norm = std::sqrt(
			std::pow(jac(0,0),2.) +	std::pow(jac(0,1),2.));
		nullspace_mat(0,0) = jac(0,0) / norm;
		nullspace_mat(1,0) = jac(0,1) / norm;
	}
	else if(dim==2) {
		predictor_update(0) = 
			jac(0,1)*jac(1,2)-jac(0,2)*jac(1,1);
		predictor_update(1) = 
			jac(0,2)*jac(1,0)-jac(0,0)*jac(1,2);
		predictor_update(2) = 
			jac(0,0)*jac(1,1)-jac(0,1)*jac(1,0);

		Vector<3,double> e1,e2,e3;
		for(int i=0; i<3; i++) {
			e1(i) = predictor_update(i);
			e2(i) = jac(0,i);
		}
		e1 /= e1.norm();
		e2 /= e2.norm();
		e3 = e1^e2;

		for(int i=0; i<3; i++) {
			nullspace_mat(i,0) = e2(i);
			nullspace_mat(i,1) = e3(i);
		}
	}

	// We renormalize the predictor_update and choose the right sign to
	// avoid a change of direction on the continuation curve.
	predictor_update /= predictor_update.norm();
	predictor_update *= (predictor_update,hc_pos-old_hc_pos)>=0 ? 1. : -1.;

	// We apply the euler step with the computed predictor_update
	old_hc_pos = hc_pos;
	hc_pos = old_hc_pos + arclength_step * predictor_update;

	// We check if we have crossed one of the target hyperplane tau=2k+1,
	// in which case we constrain the updated position to be on this plane.
	// This allows us to record a new solution
	if( std::floor((std::abs(old_hc_pos(dim)/this->L)-1.)/2.) !=
		std::floor((std::abs(hc_pos(dim)/this->L)-1.)/2.) && record_sol!=true ) {

		boundary_step = 
			( (std::round((hc_pos(dim)/this->L-1.)/2.)*2.+1)*this->L
			- old_hc_pos(dim) ) / predictor_update(dim);
		hc_pos =
			old_hc_pos + boundary_step * predictor_update;
		record_sol = true;

		if(dim==1) {
			nullspace_mat(0,0) = 1;
			nullspace_mat(1,0) = 0;
		}
		else if(dim==2) {
			nullspace_mat = 0;
			nullspace_mat(0,0) = 1;
			nullspace_mat(1,1) = 1;
		}
	}
	else 
		record_sol = false;

	// We check if we have crossed the hyperplane tau=2k,
	// in which case we constrain the updated position to be on this plane.
	// This allow to detect a possible loop
	if( std::floor(old_hc_pos(dim)/(2.*this->L)) !=
			std::floor(hc_pos(dim)/(2.*this->L)) && loop_test!=true ) {
		boundary_step = 
			( std::round(hc_pos(dim)/(2.*this->L))*2.*this->L-old_hc_pos(dim) ) /
			predictor_update(dim);
		hc_pos =
			old_hc_pos + boundary_step * predictor_update;
		loop_test = true;

		if(dim==1) {
			nullspace_mat(0,0) = 1;
			nullspace_mat(1,0) = 0;
		}
		else if(dim==2) {
			nullspace_mat = 0;
			nullspace_mat(0,0) = 1;
			nullspace_mat(1,1) = 1;
		}
	}
	else 
		loop_test = false;

	return true;
}

template <int dim>
bool HCInverseSolver<dim>::corrector_step(
		std::vector<Vector<dim,double> > &sols) {

	// We run the Newton algorithm orthogonally to the predictor update
	// to improve the current position
	bool status = true;
	int iter = 0;
	do {
		iter++;
		if(!update_reduced_jacobian_and_residual() || std::isnan(hc_pos.norm()) ||
				(iter>10 && res.norm()>newton_tol*this->L)) {
			status = false;
			break;
		}
		if(res.norm()<newton_tol*this->L) {
			break;
		}
		status = newton_update();

	} while(status);

	// If the convergence failed, we restart the predictor+corrector
	// step with a smaller arclength_step
	if(res.norm()>newton_tol*this->L) {
		arclength_step *= 0.5;
		if(record_sol) {
			arclength_step = boundary_step / 2.;
			record_sol = false;
		}
		if(loop_test) {
			arclength_step = boundary_step / 2.;
			loop_test = false;
		}
		hc_pos = old_hc_pos + arclength_step*predictor_update;

		return corrector_step(sols);
	}
	else if(status) {
		if(record_sol) {
			for(auto sol : sols) {
				if( (x-sol).norm()<1e-4*this->L) {
					return false;
				}
			}
			sols.push_back(x);
		}
		if(loop_test) {
			if( (x-x0).norm()<1e-4*this->L) {
				loop_detected = true;
				return false;
			}
		}
		if((hc_pos-old_hc_pos).norm()<1e-9)
			return false;
	}

	arclength_step = std::min(
		max_step, arclength_step*std::pow(2.,std::erf(opt_iter-iter)));

	return status;
}

template <int dim>
bool HCInverseSolver<dim>::newton_update() {

	// We compute the Newton update and apply the line search
	corrector_update = 0;
	for(int i=0; i<dim+1; i++)
		for(int j=0; j<dim; j++)
			for(int k=0; k<dim; k++)
				corrector_update(i) -=
					nullspace_mat(i,j)*jac_perp_inv(j,k)*res(k);

	double old_res_norm = res.norm();

	prev_hc_pos = hc_pos;
	hc_pos = prev_hc_pos + corrector_update;

	double step = 1;
	double ratio = 1./2.;
	int power = 0;
	while(!update_residual() && step>1e-9) {
		power++;
		step = std::pow(ratio,power);
		hc_pos = prev_hc_pos + step*corrector_update;
	}
	while( old_res_norm-res.norm() <= 0.5*step*old_res_norm && step>1e-9) {
		power++;
		step = std::pow(ratio,power);
		hc_pos = prev_hc_pos + step*corrector_update;

		if(!update_residual())
			return false;
	}

	if( std::floor((std::abs(prev_hc_pos(dim)/this->L)-1.)/2.) !=
		std::floor((std::abs(hc_pos(dim)/this->L)-1.)/2.) && record_sol!=true ) {

		hc_pos =
			prev_hc_pos +
			( (std::round((hc_pos(dim)/this->L-1.)/2.)*2.+1)*this->L
			- prev_hc_pos(dim) ) / corrector_update(dim) * corrector_update;
		record_sol = true;

		if(dim==1) {
			nullspace_mat(0,0) = 1;
			nullspace_mat(1,0) = 0;
		}
		else if(dim==2) {
			nullspace_mat = 0;
			nullspace_mat(0,0) = 1;
			nullspace_mat(1,1) = 1;
		}
	}

	return true;
}

template <int dim>
bool NHCInverseSolver<dim>::update_residual() {

	// We get the current position 
	for(int i=0; i<dim; i++)
		this->x(i) = this->hc_pos(i);

	// We get the map value and compute the homotopy residual
	if(!this->map.get_value(this->x,this->map_val)) {
		return false;
	}

	this->res =
		std::cos(PI*this->hc_pos(dim)/(2.*this->L)) * 
			(this->map_val-this->map_val_x0) +
		std::sin(PI*this->hc_pos(dim)/(2.*this->L)) *
			(this->map_val-this->y);
	return true;
}

template <int dim>
bool NHCInverseSolver<dim>::update_jacobian_and_residual() {

	// We update the residual
	if(!update_residual())
		return false;

	// We update the Jacobian 
	if(!this->map.get_gradient(this->x,this->map_grad)) {
		return false;
	}
	for(int i=0; i<dim; i++) {
		for(int j=0; j<dim; j++)
			this->jac(i,j) =
				( std::sin(PI*this->hc_pos(dim)/(2.*this->L))
				+ std::cos(PI*this->hc_pos(dim)/(2.*this->L)) ) *
					this->map_grad(j,i);
		this->jac(i,dim) = PI/(2.*this->L) * (
			- std::sin(PI*this->hc_pos(dim)/(2.*this->L)) *
				(this->map_val(i)-this->map_val_x0(i))
			+ std::cos(PI*this->hc_pos(dim)/(2.*this->L)) *
				(this->map_val(i)-this->y(i)));
	}

	return true;
}

template <int dim>
bool FPHCInverseSolver<dim>::update_residual() {

	// We get the current position 
	for(int i=0; i<dim; i++)
		this->x(i) = this->hc_pos(i);

	// We get the map value and compute the homotopy residual
	if(!this->map.get_value(this->x,this->map_val)) {
		return false;
	}

	this->res =
		std::cos(PI*this->hc_pos(dim)/(2.*this->L)) * (this->x-this->x0) +
		std::sin(PI*this->hc_pos(dim)/(2.*this->L)) * (this->map_val-this->y);
	return true;
}

template <int dim>
bool FPHCInverseSolver<dim>::update_jacobian_and_residual() {

	// We update the residual
	if(!update_residual())
		return false;

	// We update the Jacobian 
	if(!this->map.get_gradient(this->x,this->map_grad)) {
		return false;
	}
	for(int i=0; i<dim; i++) {
		for(int j=0; j<dim; j++)
			this->jac(i,j) =
				std::sin(PI*this->hc_pos(dim)/(2.*this->L)) * this->map_grad(j,i);
		this->jac(i,i) +=
				std::cos(PI*this->hc_pos(dim)/(2.*this->L));
		this->jac(i,dim) = PI/(2.*this->L) * (
			- std::sin(PI*this->hc_pos(dim)/(2.*this->L)) *
				(this->x(i)-this->x0(i))
			+ std::cos(PI*this->hc_pos(dim)/(2.*this->L)) *
				(this->map_val(i)-this->y(i)));
	}

	return true;
}

template <int dim>
bool HCInverseSolver<dim>::update_reduced_jacobian_and_residual() {

	if(!update_jacobian_and_residual())
		return false;

	// We obtain the reduced jacobian by right-multiplication with the
	// nullspace matrix
	jac_perp = 0;
	for(int i=0; i<dim; i++)
		for(int j=0; j<dim; j++)
			for(int k=0; k<dim+1; k++)
				jac_perp(i,j) +=
					jac(i,k)*nullspace_mat(k,j);

	// We inverse the reduced jacobian
	if(dim==1) {
		jac_perp_det = jac_perp(0,0);
		jac_perp_inv(0,0) = 1./jac_perp_det;
	}
	else if(dim==2) {
		jac_perp_det = jac_perp(0,0)*jac_perp(1,1)-jac_perp(1,0)*jac_perp(0,1);
		jac_perp_inv(0,0) = jac_perp(1,1) / jac_perp_det;
		jac_perp_inv(0,1) = -jac_perp(0,1) / jac_perp_det;
		jac_perp_inv(1,0) = -jac_perp(1,0) / jac_perp_det;
		jac_perp_inv(1,1) = jac_perp(0,0) / jac_perp_det;
	}

	return true;
}

template class HCInverseSolver<1>;
template class HCInverseSolver<2>;

template class NHCInverseSolver<1>;
template class NHCInverseSolver<2>;

template class FPHCInverseSolver<1>;
template class FPHCInverseSolver<2>;
