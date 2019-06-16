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
		coefs),
	wavelengths(coefs.wavelengths()) {}

void VolumeOutput::apply(
		VectorField<double> &lc_sol,
		std::vector<VectorField<std::complex<double> > > (&bpm_sol)[2]) {

	std::string filename = format_global_path() + ".vti";
	boost::filesystem::path path(filename);

	std::cout << "Exporting bulk optical field to " << path.filename().string() << std::endl;

	// First, we copy all the data in vtk structures
	unsigned int comp, i;

	unsigned int Nx = lc_sol.mesh.Nx;
	unsigned int Ny = lc_sol.mesh.Ny;
	unsigned int Nz = lc_sol.mesh.Nz;

	double delta_x = lc_sol.mesh.delta_x;
	double delta_y = lc_sol.mesh.delta_y;
	double delta_z = lc_sol.mesh.delta_z;

	auto vti_data = vtkSmartPointer<vtkImageData>::New();
	vti_data->SetDimensions(Nx, Ny, Nz);
	vti_data->SetSpacing(delta_x, delta_y, delta_z);
	vti_data->SetOrigin(-delta_x*(Nx-1)/2., -delta_y*(Ny-1)/2., -delta_z*(Nz-1)/2.);

	std::string pol_strs[2] ={"inputX_", "inputY_"};
	for(unsigned int wave_idx=0; wave_idx<wavelengths.size(); wave_idx++) {
		std::string suffix_str =
			std::to_string(wavelengths[wave_idx])+"um";

		for(unsigned int pol_idx=0; pol_idx<2; pol_idx++) {
			auto real_data = vtkSmartPointer<vtkDoubleArray>::New();
			real_data->SetNumberOfComponents(3);
			real_data->SetNumberOfTuples(Nx*Ny*Nz);
			auto varname = "E_real_"+pol_strs[pol_idx]+suffix_str;
			real_data->SetName(varname.c_str());
			
			auto imag_data = vtkSmartPointer<vtkDoubleArray>::New();
			imag_data->SetNumberOfComponents(3);
			imag_data->SetNumberOfTuples(Nx*Ny*Nz);
			varname = "E_imag_"+pol_strs[pol_idx]+suffix_str;
			imag_data->SetName(varname.c_str());

			for(i=0; i<Nx*Ny*Nz; i++) {
				for(comp=0; comp<2; comp++) {
					real_data->SetComponent(
						i, comp, std::real(bpm_sol[pol_idx][wave_idx](i,comp)));
					imag_data->SetComponent(
						i, comp, std::imag(bpm_sol[pol_idx][wave_idx](i,comp)));
				}
				real_data->SetComponent(i, 2, 0);
				imag_data->SetComponent(i, 2, 0);
			}
			vti_data->GetPointData()->AddArray(real_data);
			vti_data->GetPointData()->AddArray(imag_data);
		}
	}

	auto writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputData(vti_data);
	writer->SetDataModeToAppended();
	writer->EncodeAppendedDataOff();
	writer->Write();
}
