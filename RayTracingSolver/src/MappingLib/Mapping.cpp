#include <fstream>
#include <regex>

#include "Mapping.h"

template <int dim1, int dim2, typename T>
Mapping<dim1,dim2,T>::Mapping(
		const std::shared_ptr<DefinitionDomain<dim1> > &def_domain) :
	def_domain(def_domain) {};

template <int dim1, int dim2, typename T>
bool Mapping<dim1,dim2,T>::get_value_embedded(
		const Vector<3,double> &x, Vector<dim2,T> &res) {

	Vector<dim1,double> x_emb;
	for(int i=0; i<dim1; i++)
		x_emb(i) = x(i+3-dim1);

	return get_value(x_emb, res);
}

template <int dim1, int dim2, typename T>
bool Mapping<dim1,dim2,T>::get_gradient_embedded(
		const Vector<3,double> &x, Matrix<3,dim2,T> &res) {

	int i,j;
	Vector<dim1,double> x_emb;
	for(i=0; i<dim1; i++)
		x_emb(i) = x(i+3-dim1);

	Matrix<dim1,dim2,T> grad_emb;
	if(get_gradient(x_emb, grad_emb)) {
		res = 0;
		for(i=0; i<dim1; i++)
			for(j=0; j<dim2; j++)
				res(i+3-dim1,j) = grad_emb(i,j);
		return true;
	}
	else
		return false;
}

template <int dim1, int dim2, typename T>
InterpolatedMapping<dim1,dim2,T>::InterpolatedMapping(
		const std::shared_ptr<std::vector<Vector<dim2,T> > > &values,
		const std::shared_ptr<CartesianMesh<dim1> > &mesh,
		const std::shared_ptr<DefinitionDomain<dim1> > &def_domain,
		int pol_order) :
	Mapping<dim1,dim2,T>::Mapping(def_domain),
	values(values),
	mesh(mesh),
	pol_order(pol_order),
	last_cell_origin_idx(-1) {

	Assert(
		values.use_count()>0,
		"The pointer of the values array is empty (InterpolatedMapping)");
	Assert(
		values->size()==mesh->n_points,
		"The total number of points is different from the size "
		"of the value array (InterpolatedMapping)");

	int n = 1;
	for(int d=0; d<dim1; d++)
		n *= pol_order+1;
	pol_coefs.resize(n);

	n = 1;
	for(int d=0; d<dim1; d++)
		n *= pol_order+1;
	v.resize(n);

	get_nearby_points_indices();
	extrapolate_data();
}

template <int dim1, int dim2, typename T>
bool InterpolatedMapping<dim1,dim2,T>::get_value(
		const Vector<dim1,double> &x, Vector<dim2,T> &res) {

	int Nc = (pol_order-1)/2;

	// No computation if x is outside the definition domain
	if( !(this->def_domain->is_inside(x)) || std::isnan(x.norm()) )
		return false;

	// We compute the mesh "indices" of x
	Vector<dim1,double> x_mesh_indices =
		( x - this->mesh->origin ) / this->mesh->cell_lengths;

	// We compute the indices of the origin of the cell associated with x
	// and the renormalized coordinate of x on this cell
	MultiDimIndex<dim1> cell_origin_indices(this->mesh->n_points_per_dim);
	Vector<dim1,double> x_cell_coords;

	int d;
	for(d=0; d<dim1; ++d) {
		// We ensure that the cell associated to x is not too near the
		// mesh boundary, in which case interpolation is not possible
		if( x_mesh_indices(d)<Nc ||
				x_mesh_indices(d)>=this->mesh->n_points_per_dim(d)-1-Nc ) {
			if(x_mesh_indices(d)>=Nc-SURFACE_EPS &&
					x_mesh_indices(d)<this->mesh->n_points_per_dim(d)-1-Nc)
				x_mesh_indices(d) = Nc+SURFACE_EPS;
			else if(x_mesh_indices(d)>Nc && x_mesh_indices(d) <= 
					this->mesh->n_points_per_dim(d)-1-Nc+SURFACE_EPS)
				x_mesh_indices(d) = 
					this->mesh->n_points_per_dim(d)-1-Nc-SURFACE_EPS;
			else
				return false;
		}
		cell_origin_indices(d) = floor(x_mesh_indices(d));
	}
	x_cell_coords = x_mesh_indices - cell_origin_indices.get();
	
	bool same_cell = true;
	if(cell_origin_indices()!=last_cell_origin_idx) {
		same_cell = false;
		last_cell_origin_idx = cell_origin_indices();
		this->assemble_pol_weights(cell_origin_indices);
	}

	res = T(0);
	MultiDimIndex<dim1> pol_idx(1+pol_order);
	for(; pol_idx.valid(); pol_idx++) {
		res += pol_coefs[pol_idx()] * pow(x_cell_coords,pol_idx);
	}
	
	return true;
}

template <int dim1, int dim2, typename T>
bool InterpolatedMapping<dim1,dim2,T>::get_gradient(
		const Vector<dim1,double> &x, Matrix<dim1,dim2,T> &res) {

	int Nc = (pol_order-1)/2;

	// No computation if x is outside the definition domain
	if( !(this->def_domain->is_inside(x)) || std::isnan(x.norm()) )
		return false;

	// We compute the mesh "indices" of x
	Vector<dim1,double> x_mesh_indices =
		( x - this->mesh->origin ) / this->mesh->cell_lengths;

	// We compute the indices of the origin of the cell associated with x
	// and the renormalized coordinate of x on this cell
	MultiDimIndex<dim1> cell_origin_indices(this->mesh->n_points_per_dim);
	Vector<dim1,double> x_cell_coords;

	int d;
	for(d=0; d<dim1; ++d) {
		// We ensure that the cell associated to x is not too near the
		// mesh boundary, in which case interpolation is not possible
		if( x_mesh_indices(d)<Nc ||
				x_mesh_indices(d)>=this->mesh->n_points_per_dim(d)-1-Nc ) {
			if(x_mesh_indices(d)>=Nc-SURFACE_EPS &&
					x_mesh_indices(d)<this->mesh->n_points_per_dim(d)-1-Nc)
				x_mesh_indices(d) = Nc+SURFACE_EPS;
			else if(x_mesh_indices(d)>Nc && x_mesh_indices(d) <= 
					this->mesh->n_points_per_dim(d)-1-Nc+SURFACE_EPS)
				x_mesh_indices(d) = 
					this->mesh->n_points_per_dim(d)-1-Nc-SURFACE_EPS;
			else
				return false;
		}
		cell_origin_indices(d) = floor(x_mesh_indices(d));
	}
	x_cell_coords = x_mesh_indices - cell_origin_indices.get();

	bool same_cell = true;
	if(cell_origin_indices()!=last_cell_origin_idx) {
		same_cell = false;
		last_cell_origin_idx = cell_origin_indices();
		this->assemble_pol_weights(cell_origin_indices);
	}

	res = 0;
	MultiDimIndex<dim1> pol_idx(1+pol_order);
	MultiDimIndex<dim1> exponent(pol_idx);
	for(; pol_idx.valid(); pol_idx++) {
		for(d=0; d<dim1; d++) {
			exponent = pol_idx;
			if(exponent(d)>0) {
				exponent(d)--;
				for(int j=0; j<dim2; j++)
					res(d,j) +=
						pol_idx(d) / this->mesh->cell_lengths(d) *
						pol_coefs[pol_idx()](j) * pow(x_cell_coords,exponent);
			}
		}
	}
	
	return true;
}

template <int dim1, int dim2, typename T>
void InterpolatedMapping<dim1,dim2,T>::normalize() {
	
	double norm;
	for(auto &val : *values) {
		norm = val.norm();
		if(norm!=0)
			val = val/norm;
	}
}

template <int dim1, int dim2, typename T>
void InterpolatedMapping<dim1,dim2,T>::extrapolate_data() {

	Matrix<dim1,dim2,T> grad;
	Vector<dim1,double> pos, inside_pos, shifted_pos;
	Vector<dim2,T> val;
	Vector<dim2,T> vals[dim1];
	int signs[dim1];

	for(auto p : nearby_point_data) {
		pos = mesh->origin + p.first.get() * mesh->cell_lengths;
		inside_pos = mesh->origin + p.second.get() * mesh->cell_lengths;
		val = values->at(p.second());

		for(int d=0; d<dim1; d++) {
			shifted_pos = pos;
			shifted_pos(d) += this->mesh->cell_lengths(d);
			if(!this->def_domain->is_inside(shifted_pos)) {
				shifted_pos(d) -= 2*this->mesh->cell_lengths(d);
				if(!this->def_domain->is_inside(shifted_pos))
					vals[d] = val;
				else
					vals[d] = values->at(p.second()-this->mesh->flatten_weights(d));
				signs[d] = -1;
			}
			else {
				vals[d] = values->at(p.second()+this->mesh->flatten_weights(d));
				signs[d] = 1;
			}
		}
		for(int i=0; i<dim1; i++)
			for(int j=0; j<dim2; j++)
				grad(i,j) = 
					( vals[i](j) - val(j) ) *
					( (double) signs[i] ) / this->mesh->cell_lengths(i);

		values->at(p.first()) = val+((pos-inside_pos)*grad);
	}
}

template <int dim1, int dim2, typename T> inline
void InterpolatedMapping<dim1,dim2,T>::get_nearby_points_indices() {

	// First, we record the origin index of cells overlapping the domain
	// interface
	int Nc = (pol_order-1)/2;
	bool valid_cell_origin, interface_cell;

	std::vector<MultiDimIndex<dim1> > interface_cells_idx;
	Vector<dim1,double> pos;

	MultiDimIndex<dim1> idx(mesh->n_points_per_dim);
	for(; idx.valid(); idx++) {

		valid_cell_origin = true;
		for(int d=0; d<dim1; d++)
			if(idx.get()(d) < Nc || idx.get()(d) >= mesh->n_points_per_dim(d)-1-Nc)
				valid_cell_origin = false;
		if(!valid_cell_origin)
			continue;

		int N_outside_pts = 0;
		MultiDimIndex<dim1> cell_idx(2);
		for(; cell_idx.valid(); cell_idx++) {
			pos = mesh->origin + (idx.get()+cell_idx.get())*mesh->cell_lengths;
			if(!this->def_domain->is_inside(pos,-SURFACE_EPS))
				N_outside_pts++;
		}
		if(N_outside_pts>0 && N_outside_pts<mesh->n_points_per_cell)
			interface_cells_idx.push_back(idx);
	}

	// Then, we iterate over interface cells and find the outside points
	// contained in a pol_order*pol_order array of cells centered on each
	// interface cells
	std::vector<MultiDimIndex<dim1> > nearby_point_idx;
	for(auto idx : interface_cells_idx) {
		MultiDimIndex<dim1> cell_idx(1+pol_order, -Nc);
		for(; cell_idx.valid(); cell_idx++) {
			pos = mesh->origin + (idx.get()+cell_idx.get())*mesh->cell_lengths;
			if(!this->def_domain->is_inside(pos,SURFACE_EPS)) {
				nearby_point_idx.push_back(idx);
				nearby_point_idx.back() = idx(cell_idx);
			}
		}
	}

	// We remove duplicates with a radix sort
	std::vector<bool> bits(mesh->n_points, false);
	for(auto idx : nearby_point_idx)
		bits[idx()] = true;
	for(auto idx : nearby_point_idx) {
		if(bits[idx()]) {
			bits[idx()] = false;
			nearby_point_data.push_back(std::make_pair(idx,idx));
		}
	}

	// Finally, for each nearby point, we compute the nearest point
	// inside the definition domain (which will contain valid mapping
	// data).
	Vector<dim1,double> inside_pos, mesh_indices;
	Vector<dim1,unsigned int> cell_origin_indices;
	for(auto& p : nearby_point_data) {
		pos = mesh->origin + p.first.get()*mesh->cell_lengths;
		this->def_domain->compute_projection(pos, inside_pos, -SURFACE_EPS);

		mesh_indices =
			( inside_pos - this->mesh->origin ) / this->mesh->cell_lengths;
		for(int d=0; d<dim1; ++d) {
			cell_origin_indices(d) = floor(mesh_indices(d));
		}

		MultiDimIndex<dim1> cell_idx(2);
		double min_distance = -1;
		for(; cell_idx.valid(); cell_idx++) {
			inside_pos =
				mesh->origin +
				(cell_origin_indices+cell_idx.get())*mesh->cell_lengths;
			if(this->def_domain->is_inside(inside_pos,SURFACE_EPS) &&
					((pos-inside_pos).norm()<min_distance || min_distance<0)) {
				min_distance = (pos-inside_pos).norm();
				idx = this->mesh->flatten(cell_origin_indices+cell_idx.get());
			}
		}
		if(min_distance<0) {
			throw(std::string("Cannot get neighbor points"));
		}
		p.second = idx;
	}
}



template class Mapping<1,1,double>;
template class Mapping<1,2,double>;
template class Mapping<1,3,double>;
template class Mapping<2,1,double>;
template class Mapping<2,2,double>;
template class Mapping<2,3,double>;
template class Mapping<3,1,double>;
template class Mapping<3,2,double>;
template class Mapping<3,3,double>;

template class Mapping<1,1,std::complex<double> >;
template class Mapping<1,2,std::complex<double> >;
template class Mapping<1,3,std::complex<double> >;
template class Mapping<2,1,std::complex<double> >;
template class Mapping<2,2,std::complex<double> >;
template class Mapping<2,3,std::complex<double> >;
template class Mapping<3,1,std::complex<double> >;
template class Mapping<3,2,std::complex<double> >;
template class Mapping<3,3,std::complex<double> >;

template class InterpolatedMapping<1,1,double>;
template class InterpolatedMapping<1,2,double>;
template class InterpolatedMapping<1,3,double>;
template class InterpolatedMapping<2,1,double>;
template class InterpolatedMapping<2,2,double>;
template class InterpolatedMapping<2,3,double>;
template class InterpolatedMapping<3,1,double>;
template class InterpolatedMapping<3,2,double>;
template class InterpolatedMapping<3,3,double>;

template class InterpolatedMapping<1,1,std::complex<double> >;
template class InterpolatedMapping<1,2,std::complex<double> >;
template class InterpolatedMapping<1,3,std::complex<double> >;
template class InterpolatedMapping<2,1,std::complex<double> >;
template class InterpolatedMapping<2,2,std::complex<double> >;
template class InterpolatedMapping<2,3,std::complex<double> >;
template class InterpolatedMapping<3,1,std::complex<double> >;
template class InterpolatedMapping<3,2,std::complex<double> >;
template class InterpolatedMapping<3,3,std::complex<double> >;
