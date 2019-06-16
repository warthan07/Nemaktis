#include <complex>
#include <iostream>
#include <fstream>
#include <string>

#include <boost/regex.hpp>
#include <boost/filesystem.hpp>

#include <Eigen/Core>

#include <vtkDataArraySelection.h>
#include <vtkXMLImageDataReader.h>

#include "FullSampleSimulation.h"
#include "DropletSampleSimulation.h"

using json = nlohmann::json;

void print_usage() {
	std::cout <<
		std::endl <<
		"Usage: ./rt-solver [ -h | -c FILE | -x FILE ]" <<
		std::endl << std::endl <<
		"\t-h" << std::endl <<
		"\t\tPrint this message and exit." << 
		std::endl << std::endl <<
		"\t-c FILE" << std::endl <<
		"\t\tCreate a default settings file with the specified name and exit." << 
		std::endl << std::endl <<
		"\t-x FILE" << std::endl <<
		"\t\tRun the code with the specified settings file." << 
		std::endl << std::endl;
}

std::string minify_json(std::string filename) {

	std::ifstream f(filename);
	if(!f.is_open())
		throw std::string(
			"Could not open the specified setting file");

	boost::regex empty_re("\\s*");
	boost::regex comment_re("\\s*#.*");
	boost::smatch match_result;

	std::string line, minified_json;
	while(std::getline(f, line)) {
		if(!boost::regex_match(line, match_result, empty_re) &&
				!boost::regex_match(line, match_result, comment_re))
			minified_json.append(line+"\n");
	}
	return minified_json;
}

void create_default_settings_file(std::string filename) {

	const char *file_content = R"V0G0N(
{
	"Light source": {
		# Widths (µm) of the light source (array of size [spatial_dim-1]). Should be smaller
		# than the director field mesh in the transverse direction if the sample type is 
		# "Full".
		"Source widths": [10, 10],

		# Number of rays per dimension for the light source
		"Source N rays per dim": [200, 200],

		# Mean wavelength (µm) for the spectrum of the source
		"Mean wavelength": 0.6,

		# Full width (µm) of the light spectrum
		"Spectral FWHM": 0.2,

		# Number of wavelengths in the spectrum
		"N wavelengths": 1
	},
	"Geometry": {
		# Type of sample: "Droplet" or "Full". If "Droplet" is choosed, the given vti file for
		# the director field should have a cubic mesh (only values inside a sphere of same
		# diameter than the mesh will be considered).
		"Sample type": "Full",

		# Relative path to a VTI file containing a VTK array "n" for the director values.
		# IMPORTANT: the specified director field should have C_1 regularity (which exclude
		# the presence of defects and n->-n jumps).
		"Director field VTI file": "",

		"Droplet sample parameters": {
			# Since the ray mapping is singular near the vertical part of the droplet boundary,
			# the light source for droplet rays is shrinked by the following factor
			"Source shrink factor": 0.92,

			# Distance (µm) between the boundary of the droplet and the sample plates (the
			# droplet is centered inside the sample).
			"Distance from upper sample plate": 10
		}
	},
	"Material properties": {
		# Ordinary refractive index of the liquid crystal
		"Ordinary refractive index": 1.5,

		# Extraordinary refractive index of the liquid crystal
		"Extraordinary refractive index": 1.6,

		# Refractive index of the host fluid in case of a droplet sample
		"Host fluid refractive index": 1.55,

		# Refractive index of the external medium outside the sample
		"External medium refractive index": 1,

		# Refractive indices of the isotropic layers defining the lower part of the sample
		# (light propagates from bottom to top in the z direction, layers should be specified
		# with increasing z)
		"Refractive indices of the lower isotropic layers": [1.5],

		# Refractive indices of the isotropic layers defining the upper part of the sample
		# (light propagates from bottom to top in the z direction, layers should be specified
		# with increasing z)
		"Refractive indices of the upper isotropic layers": [1.5],

		# Thicknesses of the isotropic layers defining the lower part of the sample
		# (light propagates from bottom to top in the z direction, layers should be specified
		# with increasing z)
		"Thicknesses of the lower isotropic layers": [1000],

		# Thicknesses of the isotropic layers defining the upper part of the sample
		# (light propagates from bottom to top in the z direction, layers should be specified
		# with increasing z)
		"Thicknesses of the upper isotropic layers": [1000]
	},
	"Visualisation": {
		# Results will be stored in this folder (automatically created if it does not exists)
		"Results folder name": "results/",

		# Widths (µm) of the target plane on which fields are reconstructed. Will be used to 
		# set the transverse size of output data. Should be smaller than the source widths.
		"Target output widths": [8, 8],

		# Number of pixels per dimension for the target plane on which fields are 
		# reconstructed. Will be used to set the transverse dims of output data
		"Target N pixels per dim": [100, 100],

		"Bulk output": {
			# Base name for bulk data files
			"Base name": "bulk",

			# Should we run the homotopy continuation algorithm to reconstruct and save optical
			# fields in the bulk of the liquid crystal?
			"Export reconstructed fields": true,

			# Should we save bulk data associated with the rays?
			"Export ray data": true
		},
		"Screen output": {
			# Base name for screen data files
			"Base name": "screen",

			# Should we run the homotopy continuation algorithm to reconstruct and save optical
			# fields on the focal plane?
			"Export reconstructed fields": true,

			# Should we save data associated with the rays on the focal plane?
			"Export ray data": true,

			# Numerical aperture of the focal lens
			"Numerical aperture": 0.4
		}
	},
	"InverseScreenMap parameters": {
		# Number of adaptative refinements when inverting the ray mapping. Should be tweaked to
		# get a good compromise between accuracy and speed.
		"N refinement cycles": 4,

		# Homotopy settings used in the coarse phase of the inversion process
		"Coarse step HC parameters": {
			# Tolerance of the newton algorithm when computing corrector steps
			"Newton tolerance": 1e-8,
			
			# The homotopy step size will be adjusted so that the newton algorithm for the
			# corrector step runs in at most N iteration, where N is set in this parameter.
			"Optimum newton step number": 4,

			# Maximum size for the homotopy step size (µm)
			"Max arclength step": 0.1,

			# Maximum distance (µm) run by a homotopy path before stopping the search for
			# inverses.
			"Max arclength": 6
		},

		# Homotopy settings used in the refinement phase of the inversion process
		"Refinement step HC parameters": {
			# Tolerance of the newton algorithm when computing corrector steps
			"Newton tolerance": 1e-8,
			
			# The homotopy step size will be adjusted so that the newton algorithm for the
			# corrector step runs in at most N iteration, where N is set in this parameter.
			"Optimum newton step number": 4,

			# Maximum size for the homotopy step size (µm)
			"Max arclength step": 0.1,

			# Maximum distance (µm) run by a homotopy path before stopping the search for
			# inverses.
			"Max arclength": 6
		}
	}
})V0G0N";

	std::ofstream f(filename);
	if(!f.is_open())
		throw std::string(
			"Could not open the specified setting file, maybe you"
			"prepended a nonexistent directory?");

	f << file_content;
}

int main(int argc, char *argv[]) {

	std::string param_file_name;

	bool failed_parsing = (argc>1) ? false : true;
	int i=1;
	while(i<argc) {
		if(strcmp(argv[i],"-x")==0) {
			if(argv[++i][0]!='-') {
				param_file_name = argv[i];
				break;
			}
			else {
				failed_parsing = true;
				break;
			}
		}
		else if(strcmp(argv[i],"-h")==0) {
			print_usage();
			return 0;
		}
		else if(strcmp(argv[i],"-c")==0) {
			if(argv[++i][0]!='-') {
				create_default_settings_file(argv[i]);
				return 0;
			}
		}
		else {
			failed_parsing = true;
			break;
		}
	}
	if(failed_parsing) {
		std::cout <<
			std::endl << "Syntax error!" << std::endl;
		print_usage();
		return -1;
	}

	Eigen::initParallel();
	json j;

	try {
		j = json::parse(minify_json(param_file_name));

		std::string filename = j.at("Geometry").at("Director field VTI file");
		boost::filesystem::path path(filename);
		if(!boost::filesystem::exists(path))
			throw std::string(
				"The given vti file for the director file does not exists");

		vtkObject::GlobalWarningDisplayOff();

		std::cout <<
			"Loading director field values" << std::endl;
		auto reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
		reader->SetFileName(filename.c_str());
		reader->UpdateInformation();

		bool found_n_data = false;
		for(unsigned int i=0; i<reader->GetNumberOfPointArrays(); i++) {
			auto array_name = reader->GetPointArrayName(i);
			if(!std::strcmp(array_name, "n")) {
				found_n_data = true;
				reader->SetPointArrayStatus(array_name, 1);
			}
			else
				reader->SetPointArrayStatus(array_name, 0);
		}
		if(!found_n_data)
			throw std::string(
				"Could not find a vector field named \"n\" in the given vti file");

		reader->Update();
		auto nfield_vti_data = reader->GetOutput();
		auto nfield_vti_vals = vtkDoubleArray::FastDownCast(
			nfield_vti_data->GetPointData()->GetAbstractArray("n"));

		int dims[3];			nfield_vti_data->GetDimensions(dims);
		double spacings[3];		nfield_vti_data->GetSpacing(spacings);

		unsigned int spatial_dim;
		if(dims[2]<=1)
			throw std::string(
				"Not enough points in the z-direction");
		else if(dims[0]==1 && dims[1]==1)
			throw std::string(
				"1D samples not supported");
		else if(dims[1]==1)
			throw std::string(
				"2D samples should be specified in the YZ plane");
		else if(dims[0]==1)
			spatial_dim = 2;
		else
			spatial_dim = 3;

		double lc_thickness = (dims[2]-1)*spacings[2];
		unsigned int N_lc_steps = dims[2];

		switch(spatial_dim) {
			case 2: {
				auto n_values = std::make_shared<std::vector<Vector<3,double> > >();
				Vector<3,double> n_val;
				for(int iy=-2; iy<dims[1]+2; iy++) {
					for(int iz=-2; iz<dims[2]+2; iz++) {
						if(iy>=0 && iy<dims[1] && iz>=0 && iz<dims[2]) 
							for(int c=0; c<3; c++)
								n_val(c) = nfield_vti_vals->GetComponent(iy+dims[1]*iz, c);
						n_values->push_back(n_val);
					}
				}
				
				Vector<2,unsigned long> full_dims({dims[1]+4,dims[2]+4});
				Vector<2,double> full_origin(
					{-(dims[1]+3)*spacings[1]/2.,-(dims[2]+3)*spacings[2]/2.});
				Vector<2,double> full_lengths(
					{(dims[1]+3)*spacings[1],(dims[2]+3)*spacings[2]});
				auto full_mesh = std::make_shared<CartesianMesh<2> >(
					full_origin, full_lengths, full_dims);

				std::shared_ptr<Simulation<2> > sim;
				if(j.at("Geometry").at("Sample type") == "Full") {
					Vector<2,unsigned long> lc_dims({dims[1],dims[2]});
					Vector<2,double> lc_origin(
						{-(dims[1]-1)*spacings[1]/2.,-(dims[2]-1)*spacings[2]/2.});
					Vector<2,double> lc_lengths(
						{(dims[1]-1)*spacings[1],(dims[2]-1)*spacings[2]});
					auto lc_mesh = std::make_shared<CartesianMesh<2> >(
						lc_origin, lc_lengths, lc_dims);
					auto lc_domain = std::make_shared<ParallelotopeDomain<2> >(*lc_mesh);
					auto n_field = std::make_shared<CubicInterpolatedMapping<2,3,double> >(
						n_values, full_mesh, lc_domain);
					n_field->normalize();

					sim = std::make_shared<FullSampleSimulation<2> >(
						j, lc_thickness, N_lc_steps, n_field);
				}
				else
					throw(std::string(
						"Error: \"Sample type\" value should be \"Full\"\n"
						"in dimension 2."));
				sim->run();
				break;
			}
			case 3: {
				auto n_values = std::make_shared<std::vector<Vector<3,double> > >();
				Vector<3,double> n_val;
				for(int ix=-2; ix<dims[0]+2; ix++) {
					for(int iy=-2; iy<dims[1]+2; iy++) {
						for(int iz=-2; iz<dims[2]+2; iz++) {
							if(ix>=0 && ix<dims[0] && iy>=0 && iy<dims[1]
									&& iz>=0 && iz<dims[2]) 
								for(int c=0; c<3; c++)
									n_val(c) = nfield_vti_vals->GetComponent(
										ix+dims[0]*(iy+dims[1]*iz), c);
							else
								n_val = 0;
							n_values->push_back(n_val);
						}
					}
				}
				
				Vector<3,unsigned long> full_dims({dims[0]+4,dims[1]+4,dims[2]+4});
				Vector<3,double> full_origin({-(dims[0]+3)*spacings[0]/2.,
					-(dims[1]+3)*spacings[1]/2.,-(dims[2]+3)*spacings[2]/2.});
				Vector<3,double> full_lengths({(dims[0]+3)*spacings[0],
					(dims[1]+3)*spacings[1],(dims[2]+3)*spacings[2]});
				auto full_mesh = std::make_shared<CartesianMesh<3> >(
					full_origin, full_lengths, full_dims);

				std::shared_ptr<Simulation<3> > sim;
				if(j.at("Geometry").at("Sample type") == "Full") {
					Vector<3,unsigned long> lc_dims({dims[0],dims[1],dims[2]});
					Vector<3,double> lc_origin({-(dims[0]-1)*spacings[0]/2.,
						-(dims[1]-1)*spacings[1]/2.,-(dims[2]-1)*spacings[2]/2.});
					Vector<3,double> lc_lengths({(dims[0]-1)*spacings[0],
						(dims[1]-1)*spacings[1],(dims[2]-1)*spacings[2]});
					auto lc_mesh = std::make_shared<CartesianMesh<3> >(
						lc_origin, lc_lengths, lc_dims);
					auto lc_domain = std::make_shared<ParallelotopeDomain<3> >(*lc_mesh);
					auto n_field = std::make_shared<CubicInterpolatedMapping<3,3,double> >(
						n_values, full_mesh, lc_domain);
					n_field->normalize();


					sim = std::make_shared<FullSampleSimulation<3> >(
						j, lc_thickness, N_lc_steps, n_field);
				}
				else if(j.at("Geometry").at("Sample type") == "Droplet") {
					auto lc_domain = std::make_shared<SphericalDomain<3> >(
						Vector<3,double>({0,0,0}), (dims[2]-1)*spacings[2]/2.);
					auto n_field = std::make_shared<CubicInterpolatedMapping<3,3,double> >(
						n_values, full_mesh, lc_domain);
					n_field->extrapolate_data();
					n_field->normalize();
					lc_thickness += 2*j.at("Geometry").at("Droplet sample parameters").at(
						"Distance from upper sample plate").get<double>();

					sim = std::make_shared<DropletSampleSimulation>(
						j, lc_thickness, N_lc_steps, n_field);
				}
				else
					throw(std::string(
						"Error: \"Sample type\" value should be \"Full\"\n"
						"or \"Droplet\" in dimension 3."));
				sim->run();
				break;
			}
		}
	}
	catch(std::string &s) {
		std::cerr <<
			std::endl << std::endl <<
			"----------------------------------------------------------------" <<
			std::endl <<
			"Exception on processing: " << std::endl <<
			s << std::endl <<
			"Aborting!" << std::endl <<
			"----------------------------------------------------------------" <<
			std::endl << std::endl;
		return -1;
	}
	catch(std::exception &exc) {
		std::cerr <<
			std::endl << std::endl <<
			"----------------------------------------------------------------" <<
			std::endl <<
			"Exception on processing: " << std::endl <<
			exc.what() << std::endl <<
			"Aborting!" << std::endl <<
			"----------------------------------------------------------------" <<
			std::endl << std::endl;
		return -1;
	}
	catch(...) {
		std::cerr <<
			std::endl << std::endl <<
			"----------------------------------------------------------------" <<
			std::endl <<
			"Unknown exception!" << std::endl <<
			"Aborting!" << std::endl <<
			"----------------------------------------------------------------" <<
			std::endl << std::endl;
		return -1;
	}
	
	return 0;
}
