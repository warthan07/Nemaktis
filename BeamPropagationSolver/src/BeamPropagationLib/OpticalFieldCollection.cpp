#include "OpticalFieldCollection.h"
#include "OpticalField.h"
#include <memory>

ScreenOpticalFieldCollection::ScreenOpticalFieldCollection(
		const CartesianMesh &mesh,
		const PhysicsCoefficients &coefs) :
	mesh(mesh),
	wavelengths(coefs.wavelengths()),
	q_vals(coefs.q_vals()),
	fields(2*wavelengths.size()*q_vals.size()),
	fft_fields(2*wavelengths.size()*q_vals.size()),
	own_field_vals(true),
	stride_fields(mesh.Nx*mesh.Ny*2),
	stride_plans(stride_fields*2*q_vals.size()),
	field_vals(static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)*stride_fields*fields.size()))),
	fft_field_vals(static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)*stride_fields*fields.size()))) {

	init_fields_and_plans();
}

ScreenOpticalFieldCollection::ScreenOpticalFieldCollection(
		const CartesianMesh &mesh,
		const PhysicsCoefficients &coefs,
		std::complex<double>* user_vals,
		unsigned int n_user_vals) :
	mesh(mesh),
	wavelengths(coefs.wavelengths()),
	q_vals(coefs.q_vals()),
	fields(2*wavelengths.size()*q_vals.size()),
	fft_fields(2*wavelengths.size()*q_vals.size()),
	own_field_vals(false),
	stride_fields(mesh.Nx*mesh.Ny*2),
	stride_plans(stride_fields*2*q_vals.size()),
	field_vals(reinterpret_cast<fftw_complex*>(user_vals)),
	fft_field_vals(static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)*stride_fields*fields.size()))) {

	if(n_user_vals!=stride_fields*fields.size())
		throw std::string("Wrong dimension for the user pointer array (OpticalFieldCollection)");

	init_fields_and_plans();
}

void ScreenOpticalFieldCollection::init_fields_and_plans() {

	unsigned int iw, iq, pol_idx, flat_idx;
	for(iw=0; iw<wavelengths.size(); iw++) {
		for(iq=0; iq<q_vals.size(); iq++) {
			for(pol_idx=0; pol_idx<2; pol_idx++) {
				flat_idx = pol_idx+2*(iq+q_vals.size()*iw);
				fields[flat_idx] = std::make_shared<TransverseOpticalField>(
					mesh, reinterpret_cast<std::complex<double>*>(&field_vals[flat_idx*stride_fields]));
				fft_fields[flat_idx] = std::make_shared<TransverseOpticalField>(
					mesh, reinterpret_cast<std::complex<double>*>(&fft_field_vals[flat_idx*stride_fields]));
			}
		}
	}

	int rank = 2;
	int field_dims[] = {mesh.Ny, mesh.Nx};
	int field_size = field_dims[0]*field_dims[1];
	forward_plan = fftw_plan_many_dft(
		rank, field_dims, 4*q_vals.size(),
		field_vals, field_dims, 1, field_size,
		fft_field_vals, field_dims, 1, field_size,
		FFTW_FORWARD, FFTW_ESTIMATE);
	backward_plan = fftw_plan_many_dft(
		rank, field_dims, 4*q_vals.size(),
		fft_field_vals, field_dims, 1, field_size,
		field_vals, field_dims, 1, field_size,
		FFTW_BACKWARD, FFTW_ESTIMATE);
}


BulkOpticalFieldCollection::BulkOpticalFieldCollection(
		const CartesianMesh &mesh,
		const PhysicsCoefficients &coefs) :
	Nx(mesh.Nx),
	Ny(mesh.Ny),
	Nz(mesh.Nz),
	wavelengths(coefs.wavelengths()),
	q_vals(coefs.q_vals()) {

	double delta_x = mesh.delta_x;
	double delta_y = mesh.delta_y;
	double delta_z = mesh.delta_z;

	vti_data = vtkSmartPointer<vtkImageData>::New();
	vti_data->SetDimensions(Nx, Ny, Nz);
	vti_data->SetSpacing(delta_x, delta_y, delta_z);
	vti_data->SetOrigin(-delta_x*(Nx-1)/2., -delta_y*(Ny-1)/2., -delta_z*(Nz-1)/2.);

	auto wavelength_data = vtkSmartPointer<vtkDoubleArray>::New();
	wavelength_data->SetNumberOfComponents(1);
	wavelength_data->SetNumberOfTuples(wavelengths.size());
	wavelength_data->SetName("lambda");
	vti_data->GetFieldData()->AddArray(wavelength_data);

	auto qx_data = vtkSmartPointer<vtkDoubleArray>::New();
	qx_data->SetNumberOfComponents(1);
	qx_data->SetNumberOfTuples(q_vals.size());
	qx_data->SetName("qx");
	vti_data->GetFieldData()->AddArray(qx_data);

	auto qy_data = vtkSmartPointer<vtkDoubleArray>::New();
	qy_data->SetNumberOfComponents(1);
	qy_data->SetNumberOfTuples(q_vals.size());
	qy_data->SetName("qy");
	vti_data->GetFieldData()->AddArray(qy_data);

	std::string pol_strs[2] ={"inputX_", "inputY_"};
	for(unsigned int wave_idx=0; wave_idx<wavelengths.size(); wave_idx++) {
		wavelength_data->SetComponent(wave_idx, 0, wavelengths[wave_idx]);
		for(unsigned int q_idx=0; q_idx<q_vals.size(); q_idx++) {
			qx_data->SetComponent(q_idx, 0, q_vals[q_idx].first);
			qy_data->SetComponent(q_idx, 0, q_vals[q_idx].second);
			for(unsigned int pol_idx=0; pol_idx<2; pol_idx++) {
				real_data_arrays.push_back(vtkSmartPointer<vtkDoubleArray>::New());
				real_data_arrays.back()->SetNumberOfComponents(3);
				real_data_arrays.back()->SetNumberOfTuples(Nx*Ny*Nz);
				auto varname =
					"E_real_"+pol_strs[pol_idx]+
					std::to_string(q_idx)+"_"+std::to_string(wave_idx);
				real_data_arrays.back()->SetName(varname.c_str());

				imag_data_arrays.push_back(vtkSmartPointer<vtkDoubleArray>::New());
				imag_data_arrays.back()->SetNumberOfComponents(3);
				imag_data_arrays.back()->SetNumberOfTuples(Nx*Ny*Nz);
				varname =
					"E_imag_"+pol_strs[pol_idx]+
					std::to_string(q_idx)+"_"+std::to_string(wave_idx);
				imag_data_arrays.back()->SetName(varname.c_str());

				vti_data->GetPointData()->AddArray(real_data_arrays.back());
				vti_data->GetPointData()->AddArray(imag_data_arrays.back());
			}
		}
	}
}

void BulkOpticalFieldCollection::set_field_val(
		int wave_idx, int q_idx, int pol_idx, int comp, int mesh_idx,
		std::complex<double> val) {

	int idx = pol_idx+2*(q_idx+q_vals.size()*wave_idx);
	real_data_arrays[idx]->SetComponent(mesh_idx, comp, std::real(val));
	imag_data_arrays[idx]->SetComponent(mesh_idx, comp, std::imag(val));
}
