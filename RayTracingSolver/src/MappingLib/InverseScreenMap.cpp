#include "InverseScreenMap.h"

template <int dim>
InverseScreenMap<dim>::InverseScreenMap(
		CubicInterpolatedMapping<dim,dim,double> &_screen_map,
		const CartesianMesh<dim> &_coarse_mesh, 
		std::shared_ptr<DefinitionDomain<dim> > &_screen_domain,
		double _typical_length,
		json j) :
	screen_map(_screen_map),
	coarse_mesh(_coarse_mesh),
	fine_mesh(_coarse_mesh),
	typical_length(_typical_length),
	screen_domain(_screen_domain),
	coarse_hc_parameters(j.at("Coarse step HC parameters")), 
	refinement_hc_parameters(j.at("Refinement step HC parameters")) {

	coarse_data.resize(coarse_mesh.n_points);

	#pragma omp parallel for
	for(int i=0; i<coarse_mesh.n_points; i++) {
		MultiDimIndex<dim> idx(coarse_mesh.n_points_per_dim, 0, false);
		idx = i;

		Vector<dim,double> screen_pos = 
			coarse_mesh.origin + idx.get() * coarse_mesh.cell_lengths;

		std::vector<Vector<dim,double> > sols_1, sols_2, sols;

		try {
			if(screen_domain->is_inside(screen_pos)) {

				NHCInverseSolver<dim> nhc_solver(
						screen_map, typical_length, coarse_hc_parameters);
				nhc_solver.find_inverses(screen_pos,sols_1);

				FPHCInverseSolver<dim> fphc_solver(
						screen_map, typical_length, coarse_hc_parameters);
				fphc_solver.find_inverses(screen_pos,sols_2);
			}
		}
		catch(std::string s) {
			std::cerr <<
				std::endl << std::endl <<
				"------------------------------------------------------------" <<
				std::endl <<
				"Exception on processing: " << std::endl <<
				s << std::endl <<
				"Aborting!" << std::endl <<
				"------------------------------------------------------------" <<
				std::endl << std::endl;
		}

		sols = sols_2;
		for(auto x1 : sols_1) {
			bool record = true;
			for(auto x2 : sols_2) {
				if( (x1-x2).norm()<1e-4*typical_length) {
					record = false;
					break;
				}
			}
			if(record)
				sols.push_back(x1);
		}
		coarse_data[i] = DataPair<dim>(screen_pos,sols);
	}
}

template <int dim>
void InverseScreenMap<dim>::refine_data() {

	// We allocate the refined array and copy coarse data
	fine_mesh = coarse_mesh;
	fine_mesh.refine();
	fine_data.resize(fine_mesh.n_points);

	MultiDimIndex<dim> idx(coarse_mesh.n_points_per_dim, 0, false);
	for(; idx.valid(); idx++)
		fine_data[(fine_mesh.flatten(idx.get()*2))] = coarse_data[idx()];
	
	// If dim==2, we need to compute the refined data on the faces of the
	// coarse mesh before computing the refined data on the vertices of
	// the coarse mesh
	auto map = screen_map;
	if(dim==2) {
		#pragma omp parallel for firstprivate(map)
		for(int i=0; i<coarse_mesh.n_points; i++) {
			MultiDimIndex<dim> idx(coarse_mesh.n_points_per_dim, 0, false);
			idx = i;

			bool skip = false;
			for(int d=0; d<dim; d++) {
				if(idx.get()(d)==coarse_mesh.n_points_per_dim(d)-1) {
					skip = true;
					break;
				}
			}
			if(skip)
				continue;

			Vector<dim,int> one(1);
			int selected_idx;
			std::size_t prev_size = 0;

			MultiDimIndex<dim> cell_idx(2, 0., false);
			for(; cell_idx.valid(); cell_idx++) {
				if(coarse_data[idx(cell_idx)].second.size()>=prev_size) {
					selected_idx = idx(cell_idx);
					prev_size = coarse_data[idx(cell_idx)].second.size();
				}
			}
			
			int flat_idx = fine_mesh.flatten(2*idx.get()+one);
			fine_data[flat_idx].first =
				fine_mesh.origin + (2*idx.get()+one)*fine_mesh.cell_lengths;

			try {
				if(screen_domain->is_inside(fine_data[flat_idx].first)) {
					NHCInverseSolver<dim> nhc_solver(
						map, typical_length, refinement_hc_parameters);
					nhc_solver.find_inverses(
						fine_data[flat_idx], coarse_data[selected_idx]);
				}
				else
					fine_data[flat_idx].second =
						std::vector<Vector<dim,double> >();
			}
			catch(std::string s) {
				std::cerr <<
					std::endl << std::endl <<
					"------------------------------------------------------------" <<
					std::endl <<
					"Exception on processing: " << std::endl <<
					s << std::endl <<
					"------------------------------------------------------------" <<
					std::endl << std::endl;
			}
		}
	}

	for(int ds=0; ds<dim; ds++) {
		Vector<dim,int> shift = 0;
		shift(ds) = 1;
		
		map = screen_map;
		#pragma omp parallel for firstprivate(map)
		for(int i=0; i<coarse_mesh.n_points; i++) {
			MultiDimIndex<dim> idx(coarse_mesh.n_points_per_dim, 0, false);
			idx = i;

			Vector<dim,int> cur_idx = 2*idx.get()+shift;

			bool skip = false;
			for(int d=0; d<dim; d++) {
				if(cur_idx(d)>=fine_mesh.n_points_per_dim(d)) {
					skip = true;
					break;
				}
			}
			if(skip)
				continue;

			int flat_idx, selected_idx;
			std::size_t prev_size = 0;
			for(int d=0; d<dim; d++) {
				if(cur_idx(d)!=0) {
					flat_idx =
						fine_mesh.flatten(cur_idx) - fine_mesh.flatten_weights(d);
					if(fine_data[flat_idx].second.size()>=prev_size) {
						selected_idx = flat_idx;
						prev_size = fine_data[flat_idx].second.size();
					}
				}
				if(cur_idx(d)!=fine_mesh.n_points_per_dim(d)-1) {
					flat_idx =
						fine_mesh.flatten(cur_idx) + fine_mesh.flatten_weights(d);
					if(fine_data[flat_idx].second.size()>=prev_size) {
						selected_idx = flat_idx;
						prev_size = fine_data[flat_idx].second.size();
					}
				}
			}

			flat_idx = fine_mesh.flatten(cur_idx);
			fine_data[flat_idx].first =
				fine_mesh.origin + cur_idx*fine_mesh.cell_lengths;

			try {
				if(screen_domain->is_inside(fine_data[flat_idx].first)) {

					NHCInverseSolver<dim> nhc_solver(
						map, typical_length, refinement_hc_parameters);
					nhc_solver.find_inverses(
						fine_data[flat_idx], fine_data[selected_idx]);
				}
				else
					fine_data[flat_idx].second =
						std::vector<Vector<dim,double> >();
			}
			catch(std::string s) {
				std::cerr <<
					std::endl << std::endl <<
					"------------------------------------------------------------" <<
					std::endl <<
					"Exception on processing: " << std::endl <<
					s << std::endl <<
					"------------------------------------------------------------" <<
					std::endl << std::endl;
			}
		}
	}

	// Finally, we copy the new data on the coarse object in order to
	// be ready for a new refinement
	coarse_mesh = fine_mesh;
	coarse_data.clear();
	coarse_data = fine_data;
}

template class InverseScreenMap<1>;
template class InverseScreenMap<2>;
