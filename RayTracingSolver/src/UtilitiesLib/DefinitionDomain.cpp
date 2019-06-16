#include <cmath>

#include "DefinitionDomain.h"

template <int dim>
bool DefinitionDomain<dim>::is_inside_embedded(
		const Vector<3,double> &pos, double eps) const {

	// We compute the projection of the point on the definition
	// domain by removing the invariant directions
	Vector<dim,double> x;
	for(int i=0; i<dim; i++)
		x(i) = pos(i+3-dim);
	return is_inside(x,eps);
}

template <int dim>
bool DefinitionDomain<dim>::is_on_surface_embedded(
		const Vector<3,double> &pos, double eps) const {

	Vector<dim,double> x;
	for(int i=0; i<dim; i++)
		x(i) = pos(i+3-dim);
	return is_on_surface(x, eps);
}

template <int dim>
bool DefinitionDomain<dim>::get_normal_embedded(
		const Vector<3,double> &pos, Vector<3,double> &res, double eps) const {

	Vector<dim,double> x,nu;
	for(int i=0; i<dim; i++)
		x(i) = pos(i+3-dim);
	
	if(get_normal(x,nu,eps)) {
		res = 0;
		for(int i=0; i<dim; i++)
			res(i+3-dim) = nu(i);
		return true;
	}
	else
		return false;
}

template <int dim>
double DefinitionDomain<dim>::compute_intersection_embedded(
		const Vector<3,double> &pos, const Vector<3,double> &vel,
		Vector<3,double> &res) const {

	Vector<dim,double> x,v,r;
	for(int i=0; i<dim; i++) {
		x(i) = pos(i+3-dim);
		v(i) = vel(i+3-dim);
	}

	double s = compute_intersection(x,v,r);
	if(s>=0) {
		res = pos + vel * s;
		return s;
	}
	else
		return -1;
}

template <int dim>
double DefinitionDomain<dim>::compute_projection_embedded(
		const Vector<3,double> &pos, Vector<3,double> &res,
		double eps) const {

	Vector<dim,double> x, r;
	for(int i=0; i<dim; i++)
		x(i) = pos(i+3-dim);
	double d = compute_projection(x, r);

	res = 0;
	for(int i=0; i<dim; i++)
		res(i+3-dim) = r(i);

	return d;
}

template <int dim>
SphericalDomain<dim>::SphericalDomain(
		const Vector<dim,double> &sphere_center,
		double sphere_radius) :
	center(sphere_center),
	radius(sphere_radius) {}

template <int dim>
bool SphericalDomain<dim>::is_inside(
		const Vector<dim,double> &pos, double eps) const {

	Vector<dim,double> delta = pos-center;
	if(delta.norm()<=radius+eps)
		return true;
	else
		return false;
}

template <int dim>
bool SphericalDomain<dim>::is_on_surface(
		const Vector<dim,double> &pos, double eps) const {

	Vector<dim,double> delta = pos-center;
	if( std::abs(delta.norm()-radius)<eps ) {
		return true;
	}
	return false;
}

template <int dim>
bool SphericalDomain<dim>::get_normal(
		const Vector<dim,double> &pos, Vector<dim,double> &res, double eps) const {

	if(is_on_surface(pos,eps)) {
		Vector<dim,double> delta = pos-center;
		res = delta/delta.norm();
		return true;
	}
	else
		return false;
};

template <int dim>
double SphericalDomain<dim>::compute_intersection(
		const Vector<dim,double> &pos, const Vector<dim,double> &vel,
		Vector<dim,double> &res) const {

	Vector<dim,double> delta = pos-center;

	double delta_norm_2 = std::pow(delta.norm(), 2.);
	double vel_norm_2 = std::pow(vel.norm(), 2.);
	double delta_vel = (delta,vel);

	double discriminant =
		delta_vel*delta_vel - vel_norm_2*(delta_norm_2-radius*radius);

	if(discriminant>0) {
		double t1 = 
			( std::sqrt(discriminant) - delta_vel ) / vel_norm_2;
		double t2 = 
			(-std::sqrt(discriminant) - delta_vel ) / vel_norm_2;

		if(t1>0 && (t1<=t2 || t2<=0) ) {
			res = pos + t1*vel;
			return t1;
		}
		else if(t2>0 && (t2<=t1 || t1<=0) ) {
			res = pos + t2*vel;
			return t2;
		}
		else
			return -1;
	}
	else
		return -1;
};

template <int dim>
double SphericalDomain<dim>::compute_projection(
		const Vector<dim,double> &pos, Vector<dim,double> &res,
		double eps) const {

	res = pos-center;
	res *= (radius+eps)/res.norm();

	return (res-pos).norm();
}

HalfDiskDomain::HalfDiskDomain(
		double disk_radius) :
	radius(disk_radius) {}

bool HalfDiskDomain::is_inside(
		const Vector<2,double> &pos, double eps) const {

	if(pos.norm()<=radius+eps && pos(0)>=-eps)
		return true;
	else
		return false;
}

bool HalfDiskDomain::is_on_surface(
		const Vector<2,double> &pos, double eps) const {

	if( (std::abs(pos.norm()-radius)<eps && pos(0)>eps) ||
			(std::abs(pos(0))<eps && std::abs(pos(1))<radius+eps)) {
		return true;
	}
	return false;
}

bool HalfDiskDomain::get_normal(
		const Vector<2,double> &pos, Vector<2,double> &res, double eps) const {

	if(is_on_surface(pos,eps)) {
		if(pos(0)>0)
			res = pos/pos.norm();
		else 
			res = {-1,0};
		return true;
	}
	else
		return false;
};

double HalfDiskDomain::compute_intersection(
		const Vector<2,double> &pos, const Vector<2,double> &vel,
		Vector<2,double> &res) const {

	double delta_norm_2 = std::pow(pos.norm(), 2.);
	double vel_norm_2 = std::pow(vel.norm(), 2.);
	double delta_vel = (pos,vel);

	double discriminant =
		delta_vel*delta_vel - vel_norm_2*(delta_norm_2-radius*radius);

	if(discriminant>0) {
		double t1 = 
			( std::sqrt(discriminant) - delta_vel ) / vel_norm_2;
		double t2 = 
			(-std::sqrt(discriminant) - delta_vel ) / vel_norm_2;
		double t3 = - pos(0) / vel(0);

		if(t1>0 && (t1<=t2 || t2<=0) && 
				(t1<=t3 || t3<=0 || std::abs(pos(1)-t3*vel(1))>radius) ) {
			res = pos + t1*vel;
			return t1;
		}
		else if(t2>0 && (t2<=t1 || t1<=0) && 
				(t2<=t3 || t3<=0 || std::abs(pos(1)-t3*vel(1))>radius) ) {
			res = pos + t2*vel;
			return t2;
		}
		else if(t3>0 && (t3<=t1 || t1<=0) && (t3<=t2 || t2<=0)) {
			res = pos + t3*vel;
			return t3;
		}
		else
			return -1;
	}
	else
		return -1;
};

double HalfDiskDomain::compute_projection(
		const Vector<2,double> &pos, Vector<2,double> &res,
		double eps) const {

	if(pos(0)>0) {
		res = pos;
		res *= (radius+eps)/res.norm();
		if(pos(0)+eps<(res-pos).norm() && pos.norm()<radius+eps) 
			res = {-eps, pos(1)};
	}
	else if(std::abs(pos(1))>radius+eps) 
		res = {-eps, radius+eps};
	else 
		res = {-eps, pos(1)};

	return (res-pos).norm();
}

template <int dim>
ParallelotopeDomain<dim>::ParallelotopeDomain(
		const Vector<dim,double> &domain_origin,
		const Vector<dim,double> &domain_lengths) :
	domain_origin(domain_origin),
	domain_lengths(domain_lengths) {}

template <int dim>
ParallelotopeDomain<dim>::ParallelotopeDomain(
		const CartesianMesh<dim> &mesh) :
	domain_origin(mesh.origin),
	domain_lengths(mesh.lengths) {}

template <int dim>
bool ParallelotopeDomain<dim>::is_inside(
		const Vector<dim,double> &pos, double eps) const {

	Vector<dim,double> tmp = pos-domain_origin;
	for(int d=0; d<dim; d++) {
		if( tmp(d)<-eps || tmp(d)>=domain_lengths(d)+eps )
			return false;
	}
	return true;
}

template <int dim>
bool ParallelotopeDomain<dim>::is_on_surface(
		const Vector<dim,double> &pos, double eps) const {

	if(is_inside(pos,-eps))
		return false;

	Vector<dim,double> tmp = pos-domain_origin;
	bool on_surface = true;
	for(int d=0; d<dim; d++) {
		if(std::abs(tmp(d))>eps && std::abs(tmp(d)-domain_lengths(d))>eps &&
				(tmp(d)<=0 || tmp(d)>=domain_lengths(d)) ) {
			on_surface = false;
			break;
		}
	}
	return on_surface;
}

template <int dim>
bool ParallelotopeDomain<dim>::get_normal(
		const Vector<dim,double> &pos, Vector<dim,double> &res, double eps) const {

	Vector<dim,double> tmp = pos-domain_origin;
	if(is_on_surface(pos, eps)) {
		for(int d=0; d<dim; d++) {
			if( std::abs(tmp(d))<eps ) {
				res = 0;
				res(d) = -1;
				break;
			}
			else if( std::abs(tmp(d)-domain_lengths(d))<eps ) {
				res = 0;
				res(d) = 1;
				break;
			}
		}
		return true;
	}
	else
		return false;
};

template <int dim>
double ParallelotopeDomain<dim>::compute_intersection(
		const Vector<dim,double> &pos, const Vector<dim,double> &vel,
		Vector<dim,double> &res) const {

	Vector<dim,double> tmp = pos-domain_origin;

	double travel_time;
	double min_travel_time = -1;
	for(int d=0; d<dim; d++) {
		travel_time = 
			-tmp(d) / vel(d);
		if( travel_time>0 && is_on_surface(pos + travel_time*vel) && 
				(travel_time<min_travel_time || min_travel_time<0 ) ) {
			min_travel_time = travel_time;
			res = pos + travel_time*vel;
		}

		travel_time = 
			( domain_lengths(d) - tmp(d) ) / vel(d);
		if( travel_time>0 && is_on_surface(pos + travel_time*vel) && 
				(travel_time<min_travel_time || min_travel_time<0 ) ) {
			min_travel_time = travel_time;
			res = pos + travel_time*vel;
		}
	}
	return min_travel_time;
};

template <int dim>
double ParallelotopeDomain<dim>::compute_projection(
		const Vector<dim,double> &pos, Vector<dim,double> &res,
		double eps) const {

	Vector<dim,double> tmp = pos-domain_origin;

	double distance, min_distance = -1;

	for(int d=0; d<dim; d++) {
		Vector<dim,double> v;
		v(d) = 1;

		distance = -tmp(d);
		if( distance>0 && is_on_surface(pos + distance*v) && 
				(distance<min_distance || min_distance<0 ) ) {
			min_distance = distance;
			res = pos + (distance-eps)*v;
		}

		distance = domain_lengths(d) - tmp(d);
		if( distance>0 && is_on_surface(pos + distance*v) && 
				(distance<min_distance || min_distance<0 ) ) {
			min_distance = distance;
			res = pos + (distance+eps)*v;
		}
	}

	if(min_distance<0) {
		for(int d=0; d<dim; d++) {
			if(tmp(d)<0)
				res(d) = domain_origin(d)-eps;
			else if(tmp(d)-domain_lengths(d)>0)
				res(d) = domain_origin(d)+domain_lengths(d)+eps;
			else
				res(d) = domain_origin(d)+tmp(d);
		}
		min_distance = (tmp-res).norm();
	}

	return min_distance;
}

template <int dim>
SlabDomain<dim>::SlabDomain(
		const Vector<dim,double> &slab_normal,	double h1, double h2) :
	slab_normal(slab_normal) {

	if(h1<=h2) {
		h_1 = h1;
		h_2 = h2;
	}
	else {
		h_1 = h2;
		h_2 = h1;
	}
}

template <int dim>
bool SlabDomain<dim>::is_inside(
		const Vector<dim,double> &pos, double eps) const {

	double pos_nu = (pos,slab_normal);
	if( pos_nu<=h_2+eps && pos_nu>=h_1-eps )
		return true;
	else
		return false;
}

template <int dim>
bool SlabDomain<dim>::is_on_surface(
		const Vector<dim,double> &pos, double eps) const {

	double pos_nu = (pos,slab_normal);
	if( std::abs(pos_nu-h_2)<eps || std::abs(pos_nu-h_1)<eps )
		return true;
	else
		return false;
}

template <int dim>
bool SlabDomain<dim>::get_normal(
		const Vector<dim,double> &pos, Vector<dim,double> &res, double eps) const {

	double pos_nu = (pos,slab_normal);
	if( std::abs(pos_nu-h_2)<eps ) {
		res = slab_normal;
		return true;
	}
	else if( std::abs(pos_nu-h_1)<eps ) {
		res = -slab_normal;
		return true;
	}
	else
		return false;
}

template <int dim>
double SlabDomain<dim>::compute_intersection(
		const Vector<dim,double> &pos, const Vector<dim,double> &vel,
		Vector<dim,double> &res) const {

	double pos_nu = (pos,slab_normal);
	double vel_nu = (vel,slab_normal);
	double travel_time;

	if( ( vel_nu>0 && pos_nu<h_1 ) || 
		( vel_nu<0 && pos_nu>h_1 && pos_nu<=h_2 ) ) {

		travel_time = (h_1 - pos_nu) / vel_nu;
		res = pos + travel_time*vel;
		return travel_time;
	}
	else if( ( vel_nu<0 && pos_nu>h_2 ) || 
		( vel_nu>0 && pos_nu>=h_1 && pos_nu<h_2 ) ) {

		travel_time = (h_2 - pos_nu) / vel_nu;
		res = pos + travel_time*vel;
		return travel_time;
	}
	else
		return -1;
}

template <int dim>
double SlabDomain<dim>::compute_projection(
		const Vector<dim,double> &pos, Vector<dim,double> &res,
		double eps) const {

	double pos_nu = (pos,slab_normal);

	if(std::abs(pos_nu-h_1)<std::abs(pos_nu-h_2)) {
		res = pos + (h_1-eps-pos_nu)*slab_normal;
		return std::abs(pos_nu-h_1);
	}
	else {
		res = pos + (h_2+eps-pos_nu)*slab_normal;
		return std::abs(pos_nu-h_2);
	} 
}

template <int dim>
SubstractedDomain<dim>::SubstractedDomain(
		const std::shared_ptr<DefinitionDomain<dim> > &outer_domain,
		const std::shared_ptr<DefinitionDomain<dim> > &inner_domain) :
	outer_domain(outer_domain),
	inner_domain(inner_domain) {}

template <int dim>
bool SubstractedDomain<dim>::is_inside(
		const Vector<dim,double> &pos, double eps) const {

	if(outer_domain->is_inside(pos,eps) && !inner_domain->is_inside(pos,-eps))
		return true;
	else
		return false;
}

template <int dim>
bool SubstractedDomain<dim>::is_on_surface(
		const Vector<dim,double> &pos, double eps) const {

	if(outer_domain->is_on_surface(pos,eps) ||
			inner_domain->is_on_surface(pos,eps))
		return true;
	else 
		return false;
}

template <int dim>
bool SubstractedDomain<dim>::get_normal(
		const Vector<dim,double> &pos, Vector<dim,double> &res, double eps) const {

	if(outer_domain->get_normal(pos,res,eps))
		return true;
	else if(inner_domain->get_normal(pos,res,eps)) {
		res = -res;
		return true;
	}
	else 
		return false;
}

template <int dim>
double SubstractedDomain<dim>::compute_intersection(
		const Vector<dim,double> &pos, const Vector<dim,double> &vel,
		Vector<dim,double> &res) const {

	double s1, s2;
	Vector<dim,double> res1, res2;

	s1 = outer_domain->compute_intersection(pos,vel,res1);
	s2 = inner_domain->compute_intersection(pos,vel,res2);

	if(s1<0 && s2<0)
		return -1;
	else if(s2>0 && (s1<0 || s2<=s1)) {
		res = res2;
		return s2;
	}
	else {
		res = res1;
		return s1;
	}
}

template <int dim>
double SubstractedDomain<dim>::compute_projection(
		const Vector<dim,double> &pos, Vector<dim,double> &res,
		double eps) const {

	double d1, d2;
	Vector<dim,double> res1, res2;

	d1 = outer_domain->compute_projection(pos,res1,eps);
	d2 = inner_domain->compute_projection(pos,res2,-eps);

	if(d2<d1) {
		res = res2;
		return d2;
	}
	else {
		res = res1;
		return d1;
	}
}

template class DefinitionDomain<1>;
template class DefinitionDomain<2>;
template class DefinitionDomain<3>;

template class SphericalDomain<1>;
template class SphericalDomain<2>;
template class SphericalDomain<3>;

template class ParallelotopeDomain<1>;
template class ParallelotopeDomain<2>;
template class ParallelotopeDomain<3>;

template class SlabDomain<1>;
template class SlabDomain<2>;
template class SlabDomain<3>;

template class SubstractedDomain<1>;
template class SubstractedDomain<2>;
template class SubstractedDomain<3>;
