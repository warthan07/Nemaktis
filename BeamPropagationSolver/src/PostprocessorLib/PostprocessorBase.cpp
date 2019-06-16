#include <boost/regex.hpp>
#include <boost/filesystem.hpp>

#include <vtkDoubleArray.h>
#include <vtkPointData.h>

#include "PostprocessorBase.h"
#include "error.h"

PostprocessorBase::PostprocessorBase(
		const BasePostprocessorSettings &pp_settings,
		const std::string &results_dirname,
		const PhysicsCoefficients &coefs) :
	base_name(pp_settings.base_name),
	results_dirname(results_dirname),
	coefs(coefs) {

	Assert(
		pp_settings.activate,
		"You are trying to create a postprocessor which should not be activated");
    if(results_dirname!="" && results_dirname!=".")
    	boost::filesystem::create_directory(results_dirname);
}

std::string PostprocessorBase::format_global_path(std::string suffix) {

	std::ostringstream fstr;

    if(results_dirname!="" && results_dirname!=".") {
	    fstr << results_dirname << "/";
	    boost::filesystem::create_directory(fstr.str());
    }
	fstr << base_name << suffix;
	return fstr.str();
}
