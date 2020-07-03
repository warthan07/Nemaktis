#include <boost/filesystem.hpp>

#include <vtkXMLImageDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkPointData.h>

#include "VolumeOutput.h"

VolumeOutput::VolumeOutput(
		const RootSettings &settings,
		const PhysicsCoefficients &coefs) :
	PostprocessorBase::PostprocessorBase(
		settings.postprocessor.volume_output,
		settings.algorithm.general.results_folder_name,
		coefs) {}

void VolumeOutput::apply(BulkOpticalFieldCollection &bulk_optical_fields) {

	std::string filename = format_global_path() + ".vti";
	boost::filesystem::path path(filename);

	std::cout << "Exporting bulk optical field to " << path.filename().string() << std::endl;

	auto writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetDataModeToAppended();
	writer->EncodeAppendedDataOff();
	writer->SetInputData(bulk_optical_fields.get_vti_data());
	writer->Write();
}
