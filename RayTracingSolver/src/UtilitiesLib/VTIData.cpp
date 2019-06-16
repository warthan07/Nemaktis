#include <boost/filesystem.hpp>
#include <vtkXMLImageDataWriter.h>

#include "VTIData.h"

RayData::RayData(
		unsigned long Nx, unsigned long Ny, unsigned long Nz,
		double dx, double dy, double dz) {

	vti_data = vtkSmartPointer<vtkImageData>::New();
	vti_data->SetDimensions(Nx, Ny, Nz);
	vti_data->SetSpacing(dx, dy, dz);
	vti_data->SetOrigin(
		-dx*(Nx-1)/2., -dy*(Ny-1)/2., -dz*(Nz-1)/2.);

	deflection = vtkSmartPointer<vtkDoubleArray>::New();
	deflection->SetNumberOfComponents(2);
	deflection->SetNumberOfTuples(Nx*Ny*Nz);
	deflection->SetName("Ray deflection");
	vti_data->GetPointData()->AddArray(deflection);

	opt_length = vtkSmartPointer<vtkDoubleArray>::New();
	opt_length->SetNumberOfComponents(1);
	opt_length->SetNumberOfTuples(Nx*Ny*Nz);
	opt_length->SetName("Optical length");
	vti_data->GetPointData()->AddArray(opt_length);

	moment = vtkSmartPointer<vtkDoubleArray>::New();
	moment->SetNumberOfComponents(3);
	moment->SetNumberOfTuples(Nx*Ny*Nz);
	moment->SetName("Renormalized wavevector");
	vti_data->GetPointData()->AddArray(moment);

	std::string pol_strs[2] ={"_inputX", "_inputY"};
	for(int pol_idx=0; pol_idx<2; pol_idx++) {
		ampl[pol_idx] = vtkSmartPointer<vtkDoubleArray>::New();
		ampl[pol_idx]->SetNumberOfComponents(3);
		ampl[pol_idx]->SetNumberOfTuples(Nx*Ny*Nz);
		auto name = "Ampl"+pol_strs[pol_idx];
		ampl[pol_idx]->SetName(name.c_str());
		vti_data->GetPointData()->AddArray(ampl[pol_idx]);
	}
}

void RayData::write(std::string basedir, std::string filename) {

	boost::filesystem::path path(basedir);
	if(!boost::filesystem::exists(path))
		boost::filesystem::create_directories(path);

	auto writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
	writer->SetFileName((basedir+"/"+filename).c_str());
	writer->SetInputData(vti_data);
	writer->SetDataModeToAppended();
	writer->EncodeAppendedDataOff();
	writer->Write();
}

FieldsData::FieldsData(
		unsigned long Nx, unsigned long Ny, unsigned long Nz,
		double dx, double dy, double dz,
		std::vector<double> &wavelengths) {

	vti_data = vtkSmartPointer<vtkImageData>::New();
	vti_data->SetDimensions(Nx, Ny, Nz);
	vti_data->SetSpacing(dx, dy, dz);
	vti_data->SetOrigin(
		-dx*(Nx-1)/2., -dy*(Ny-1)/2., -dz*(Nz-1)/2.);

	std::string pol_strs[2] ={"inputX_", "inputY_"};
	for(auto wavelength : wavelengths) {
		for(int pol_idx=0; pol_idx<2; pol_idx++) {
			E_real[pol_idx].push_back(vtkSmartPointer<vtkDoubleArray>::New());
			E_real[pol_idx].back()->SetNumberOfComponents(3);
			E_real[pol_idx].back()->SetNumberOfTuples(Nx*Ny*Nz);
			auto name = "E_real_"+pol_strs[pol_idx]+std::to_string(wavelength)+"um";
			E_real[pol_idx].back()->SetName(name.c_str());
			vti_data->GetPointData()->AddArray(E_real[pol_idx].back());

			E_imag[pol_idx].push_back(vtkSmartPointer<vtkDoubleArray>::New());
			E_imag[pol_idx].back()->SetNumberOfComponents(3);
			E_imag[pol_idx].back()->SetNumberOfTuples(Nx*Ny*Nz);
			name = "E_imag_"+pol_strs[pol_idx]+std::to_string(wavelength)+"um";
			E_imag[pol_idx].back()->SetName(name.c_str());
			vti_data->GetPointData()->AddArray(E_imag[pol_idx].back());
		}
	}
	ray_multiplicity = vtkSmartPointer<vtkDoubleArray>::New();
	ray_multiplicity->SetNumberOfComponents(1);
	ray_multiplicity->SetNumberOfTuples(Nx*Ny*Nz);
	ray_multiplicity->SetName("Ray multiplicity");
	vti_data->GetPointData()->AddArray(ray_multiplicity);
}

void FieldsData::write(std::string basedir, std::string filename) {

	boost::filesystem::path path(basedir);
	if(!boost::filesystem::exists(path))
		boost::filesystem::create_directories(path);

	auto writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
	writer->SetFileName((basedir+"/"+filename).c_str());
	writer->SetInputData(vti_data);
	writer->SetDataModeToAppended();
	writer->EncodeAppendedDataOff();
	writer->Write();
}
