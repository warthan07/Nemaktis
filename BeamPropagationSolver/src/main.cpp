#include <iostream>
#include <fstream>
#include <signal.h>
#include <boost/regex.hpp>

#include "RootSettings.h"

void print_usage() {
	std::cout <<
		std::endl <<
		"Usage: ./bpm-solver [ -h | -c FILE | -x FILE ]" <<
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

void create_default_settings_file(std::string filename) {

	const char* file_content = R"V0G0N({
	"Algorithm settings": {
		"General": {
			# The type of field for the LC field: "Director" (vector of dimension 3) or
			# "Q-tensor" (vector of dimension 6 containing the {xx,yy,zz,xy,xz,yz} components,
			# since the Q-tensor is symmetric).
			"LC field type":		"Director",

    		# The name of the folder where the results will be stored.
    		# The folder will be created if it does not exists.
			"Results folder name":	"results/"
		},
		"Beam propagation": {
			# Number of Woodbury iterations (for most systems, 1 or 2
			# iterations should be enough to get accurate results, but
			# you should always try higher iterations just to check the
			# validity of the convergence)
			"N Woodbury steps": 		2,

			# Number of evolution substeps per z-slab
			"Number of substeps per slab": 1
		}
	},

	"Physics settings": {
		"Initial conditions": {
			# The type of beam profile (GaussianBeam(WAIST) or UniformBeam)
			# used to initialize the optical fields on the entrance plane.
			"Beam profile":		"UniformBeam",

			# Relative path to a vti file containing the LC orientational field, named "n"
			# in "Director" mode or "Q" in "Q-tensor" mode.
			"LC field file":	"",

			# Basis convention for the input vti file (in our code, Z is
			# the direction where light is propagated): XYZ, YZX or ZXY
			"Basis convention":	"XYZ"
		},

		"Coefficients": {
    		# Ordinary refractive index (can include a dependence on the
			# wavelength "lambda" in µm)
			"no": 				"1.5",

    		# Extraordinary refractive index (can include a dependence on the
			# wavelength "lambda" in µm)
			"ne": 				"1.6",

    		# Refractive index of the host fluid, used in case of a
			# non-trivial LC domain mask (can include a dependence on the
			# wavelength "lambda" in µm). Only used in backend-mode.
			"nhost": 				"1.5",

    		# Refractive index of the input medium below the LC layer (can
			# include a dependence on the wavelength "lambda" in µm). 
			"nin": 				"1.5",

			# Mean wavelength for the spectrum of the input beam
			"Mean wavelength":	0.6,

			# Difference between the maximum and minimum wavelength in the
			# spectrum of the input beam
			"Spectral FWHM":	0.2,

			# Number of wavelengths in the spectrum of the input beam
			"N wavelengths":	1
		}
	},

	"Postprocessor settings": {
		"Bulk output": {
    		# Do we need to run this postprocessor?
    		"Activate":								false,

    		# Base name for the file(s) exported by this postprocessor
    		"Base name":							"bulk"
		},

		"Screen output": {
    		# Do we need to run this postprocessor?
    		"Activate":								true,

    		# Base name for the file(s) exported by this postprocessor
    		"Base name":							"screen",

			# Array containing the thicknesses of all isotropic layers
			# between the LC layer and the microscope objective, in the
			# same order as they are encountered by the light beam.
			"Isotropic layer thicknesses": 			[1000],

			# Array containing the refractive indices of all isotropic layers
			# between the LC layer and the microscope objective, in the
			# same order as they are encountered by the light beam.
			"Isotropic layer refractive indices":	[1.5],

			# Vertical position of the focalisation plane with respect
			# to the central focalisation plane associated with the
			# output plane of the LC layer.
			"Focalisation z-shift":					0,

			# Numerical aperture of the focal lens
			"Numerical aperture":					0.3
		}
	}
})V0G0N";

	std::ofstream f(filename);
	if(!f.is_open())
		throw std::string(
			"Could not open the specified setting file, maybe you prepended "
			"a non-existent directory?");
	
	f << file_content;
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

int main(int argc, char *argv[]) {

	std::string param_file_name;
	std::string log_file_name = "log.txt";

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
		return 0;
	}

	try {
		auto j = nlohmann::json::parse(minify_json(param_file_name));
		RootSettings settings(j);

		std::string result_folder =
			settings.algorithm.general.results_folder_name;
		std::string create_dir = "mkdir -p " + result_folder;
		system(create_dir.c_str());

		runFromSettings(settings);
	}
	catch (std::exception &exc) {
		std::cerr <<
			std::endl << std::endl <<
			"-------------------------------------------------" <<
			std::endl <<
			"Exception on processing: " << std::endl <<
			exc.what() << std::endl <<
			"Aborting!" << std::endl <<
			"-------------------------------------------------" <<
			std::endl;
		return 1;
	}
	catch(std::string &str) {
		std::cerr <<
			std::endl << std::endl <<
			"-------------------------------------------------" <<
			std::endl <<
			"Exception on processing: " << std::endl <<
			str << std::endl <<
			"Aborting!" << std::endl <<
			"-------------------------------------------------" <<
			std::endl;
		return 1;
	}
	catch (...) {
		std::cerr <<
			std::endl << std::endl <<
			"-------------------------------------------------" <<
			std::endl <<
			"Unknown exception!" << std::endl <<
			"Aborting!" << std::endl <<
			"-------------------------------------------------" <<
			std::endl;
		return 1;
	}
	return 0;
}
