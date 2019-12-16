#include <boost/filesystem.hpp>

#include <fftw3.h>

#include <vtkXMLImageDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkPointData.h>

#include "ScreenOutput.h"

ScreenOutput::ScreenOutput(
		const RootSettings &settings,
		const PhysicsCoefficients &coefs) :
	PostprocessorBase::PostprocessorBase(
		settings.postprocessor.micrograph_output,
		settings.algorithm.general.results_folder_name,
		coefs),
	iso_layer_thickness(settings.postprocessor.micrograph_output.iso_layer_thickness),
	iso_layer_index(settings.postprocessor.micrograph_output.iso_layer_index),
	focalisation_z_shift(settings.postprocessor.micrograph_output.focalisation_z_shift),
	numerical_aperture(settings.postprocessor.micrograph_output.numerical_aperture),
	wavelengths(coefs.wavelengths()) {

	if(iso_layer_thickness.size()!=iso_layer_index.size())
		throw std::string(
			"Inconsistent number of isotropic layers");

	// index of the air between the last iso layer and the objective
	iso_layer_index.push_back(1);
}

void ScreenOutput::apply(
		VectorField<double> &lc_sol,
		std::vector<VectorField<std::complex<double> > > (&bpm_sol)[2]) {

	std::cout << "Calculating output transverse fields..." << std::endl;

	Nx = lc_sol.mesh.Nx;
	Ny = lc_sol.mesh.Ny;
	Nz = lc_sol.mesh.Nz;

	delta_x = lc_sol.mesh.delta_x;
	delta_y = lc_sol.mesh.delta_y;

	// We create the FFTW plans for the forward and backward transform
	std::cout << "    Accumulating FFTW wisdom..." << std::endl;

	fftw_complex *in = static_cast<fftw_complex*>(
		fftw_malloc(sizeof(fftw_complex) * Nx*Ny));
	fftw_complex *out = static_cast<fftw_complex*>(
		fftw_malloc(sizeof(fftw_complex) * Nx*Ny));
	fftw_plan forward_plan =
		fftw_plan_dft_2d(Ny, Nx, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan backward_plan =
		fftw_plan_dft_2d(Ny, Nx, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

	// Array of pointers for the vti data structures
	std::vector<vtkSmartPointer<vtkDoubleArray> > vti_data_arrays(4*wavelengths.size());

	std::cout <<
		"    Filtering optical fields..." << std::endl;
	std::string pol_strs[2] ={"inputX_", "inputY_"};

	#pragma omp parallel for
	for(int wave_idx=0; wave_idx<wavelengths.size(); wave_idx++) {
		// FFTW arrays for this thread
		fftw_complex *field = static_cast<fftw_complex*>(
			fftw_malloc(sizeof(fftw_complex) * Nx*Ny));
		fftw_complex *fft_field = static_cast<fftw_complex*>(
			fftw_malloc(sizeof(fftw_complex) * Nx*Ny));

		// Propagation filter for this wavelength
		auto iso_filter = assemble_iso_filter(wavelengths[wave_idx]);

		// VTK arrays for the current wavelength
		std::string suffix_str =
			std::to_string(wavelengths[wave_idx])+"um";

		for(int pol_idx=0; pol_idx<2; pol_idx++) {
			auto real_data = vtkSmartPointer<vtkDoubleArray>::New();
			real_data->SetNumberOfComponents(3);
			real_data->SetNumberOfTuples(Nx*Ny);
			auto varname = "E_real_"+pol_strs[pol_idx]+suffix_str;
			real_data->SetName(varname.c_str());
			
			auto imag_data = vtkSmartPointer<vtkDoubleArray>::New();
			imag_data->SetNumberOfComponents(3);
			imag_data->SetNumberOfTuples(Nx*Ny);
			varname = "E_imag_"+pol_strs[pol_idx]+suffix_str;
			imag_data->SetName(varname.c_str());

			vti_data_arrays[4*wave_idx+2*pol_idx] = real_data;
			vti_data_arrays[4*wave_idx+2*pol_idx+1] = imag_data;

			for(int comp=0; comp<2; comp++) {
				// First, we fill the fftw input array with our data
				for(int i=0; i<Nx*Ny; i++) {
					field[i][0] = 
						std::real(bpm_sol[pol_idx][wave_idx](i+Nx*Ny*(Nz-1),comp));
					field[i][1] = 
						std::imag(bpm_sol[pol_idx][wave_idx](i+Nx*Ny*(Nz-1),comp));
				}

				// We switch to fourrier space and apply the iso filter
				fftw_execute_dft(forward_plan, field, fft_field);
				for(int i=0; i<Nx*Ny; i++) {
					auto val = std::complex<double>(fft_field[i][0], fft_field[i][1]);
					fft_field[i][0] = std::real(val * iso_filter->at(i));
					fft_field[i][1] = std::imag(val * iso_filter->at(i));
				}

				// We come back to real space and save the transformed
				// data
				fftw_execute_dft(backward_plan, fft_field, field);

				for(int i=0; i<Nx*Ny; i++) {
					real_data->SetComponent(i, comp, field[i][0]/(Nx*Ny));
					imag_data->SetComponent(i, comp, field[i][1]/(Nx*Ny));
				}
			}
			for(int i=0; i<Nx*Ny; i++) {
				real_data->SetComponent(i, 2, 0);
				imag_data->SetComponent(i, 2, 0);
			}
		}

		// We free the FFTW pointers
		fftw_free(field);
		fftw_free(fft_field);
	}

	// We free the FFTW pointers
	fftw_destroy_plan(forward_plan);
	fftw_destroy_plan(backward_plan);
	fftw_free(in);
	fftw_free(out);

	// We export all the calculated data
	auto vti_data = vtkSmartPointer<vtkImageData>::New();
	vti_data->SetDimensions(Nx, Ny, 1);
	vti_data->SetSpacing(delta_x, delta_y, 0);
	vti_data->SetOrigin(-delta_x*(Nx-1)/2., -delta_y*(Ny-1)/2., 0);

	for(auto v : vti_data_arrays)
		vti_data->GetPointData()->AddArray(v);

	std::string filename = format_global_path() + ".vti";
	boost::filesystem::path path(filename);

	std::cout << "    Exporting " << path.filename().string() << std::endl;

	auto writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputData(vti_data);
	writer->SetDataModeToAppended();
	writer->EncodeAppendedDataOff();
	writer->Write();
}

void ScreenOutput::apply_no_export(
		VectorField<double> &lc_sol,
		std::vector<VectorField<std::complex<double> > > (&bpm_sol)[2],
		std::complex<double>* fields_vals) {

	std::cout << "Calculating output transverse fields..." << std::endl;

	Nx = lc_sol.mesh.Nx;
	Ny = lc_sol.mesh.Ny;
	Nz = lc_sol.mesh.Nz;

	delta_x = lc_sol.mesh.delta_x;
	delta_y = lc_sol.mesh.delta_y;

	// We create the FFTW plans for the forward and backward transform
	std::cout << "    Accumulating FFTW wisdom..." << std::endl;

	fftw_complex *in = static_cast<fftw_complex*>(
		fftw_malloc(sizeof(fftw_complex) * Nx*Ny));
	fftw_complex *out = static_cast<fftw_complex*>(
		fftw_malloc(sizeof(fftw_complex) * Nx*Ny));
	fftw_plan forward_plan =
		fftw_plan_dft_2d(Ny, Nx, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan backward_plan =
		fftw_plan_dft_2d(Ny, Nx, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

	std::cout <<
		"    Filtering optical fields..." << std::endl;

	#pragma omp parallel for
	for(int wave_idx=0; wave_idx<wavelengths.size(); wave_idx++) {
		// FFTW arrays for this thread
		fftw_complex *field = static_cast<fftw_complex*>(
			fftw_malloc(sizeof(fftw_complex) * Nx*Ny));
		fftw_complex *fft_field = static_cast<fftw_complex*>(
			fftw_malloc(sizeof(fftw_complex) * Nx*Ny));

		// Propagation filter for this wavelength
		auto iso_filter = assemble_iso_filter(wavelengths[wave_idx]);

		for(int pol_idx=0; pol_idx<2; pol_idx++) {
			for(int comp=0; comp<2; comp++) {
				// First, we fill the fftw input array with our data
				for(int i=0; i<Nx*Ny; i++) {
					field[i][0] = 
						std::real(bpm_sol[pol_idx][wave_idx](i+Nx*Ny*(Nz-1),comp));
					field[i][1] = 
						std::imag(bpm_sol[pol_idx][wave_idx](i+Nx*Ny*(Nz-1),comp));
				}

				// We switch to fourrier space and apply the iso filter
				fftw_execute_dft(forward_plan, field, fft_field);
				for(int i=0; i<Nx*Ny; i++) {
					auto val = std::complex<double>(fft_field[i][0], fft_field[i][1]);
					fft_field[i][0] = std::real(val * iso_filter->at(i));
					fft_field[i][1] = std::imag(val * iso_filter->at(i));
				}

				// We come back to real space and save the transformed
				// data
				fftw_execute_dft(backward_plan, fft_field, field);

				for(int i=0; i<Nx*Ny; i++)
					fields_vals[i+Nx*Ny*(comp+2*pol_idx+4*wave_idx)] =
						std::complex<double>(field[i][0]/(Nx*Ny), field[i][1]/(Nx*Ny));
			}
		}

		// We free the FFTW pointers
		fftw_free(field);
		fftw_free(fft_field);
	}

	// We free the FFTW pointers
	fftw_destroy_plan(forward_plan);
	fftw_destroy_plan(backward_plan);
	fftw_free(in);
	fftw_free(out);
}

std::shared_ptr<std::vector<std::complex<double> > > ScreenOutput::assemble_iso_filter(
		double wavelength) const {

	double z_foc = focalisation_z_shift;
	double z_end = 0;
	for(int i=0; i<iso_layer_thickness.size(); i++) {
		z_foc += iso_layer_thickness[i] * (1-1./iso_layer_index[i]);
		z_end += iso_layer_thickness[i];
	}

	auto iso_filter = std::make_shared<std::vector<std::complex<double> > >(Nx*Ny);
	double k0 = 2*PI/wavelength;

	// transfer matrix in column-major order
	std::vector<std::complex<double> > tmat(4), prev_tmat(4); 
	std::complex<double> t1, t2;

	#pragma omp parallel for firstprivate(tmat, prev_tmat, t1, t2)
	for(int iy=0; iy<Ny; iy++) {
		for(int ix=0; ix<Nx; ix++) {
			double kx = (ix<double(Nx/2.)) ?
				2*PI*ix/(delta_x*Nx) :
				-2*PI*(Nx-ix)/(delta_x*Nx);
			double ky = (iy<double(Ny/2.)) ?
				2*PI*iy/(delta_y*Ny) :
				-2*PI*(Ny-iy)/(delta_y*Ny);
			double k = std::sqrt(kx*kx+ky*ky);

			if(k<k0*numerical_aperture) {
				tmat[0] = 1;
				tmat[1] = 0;
				tmat[2] = 0;
				tmat[3] = 1;

				prev_tmat = tmat;

				double kN = std::sqrt(k0*k0-k*k);

				for(int p=0; p<iso_layer_thickness.size(); p++) {
					double np = iso_layer_index[p];
					// double np = 1.5046+0.00420/std::pow(wavelength,2.);
					double ep = iso_layer_thickness[p];
					double kp = std::sqrt(std::pow(k0*np,2.)-k*k);
					double next_kp = std::sqrt(std::pow(k0*iso_layer_index[p+1],2.)-k*k);

					t1 = 0.5*(1+kp/next_kp)*std::exp(std::complex<double>(0,(kp-kN)*ep));
					t2 = 0.5*(1-kp/next_kp)*std::exp(std::complex<double>(0,(kp-kN)*ep));

					tmat[0] = t1*prev_tmat[0] + std::conj(t2)*prev_tmat[1];
					tmat[1] = t2*prev_tmat[0] + std::conj(t1)*prev_tmat[1];
					tmat[2] = t1*prev_tmat[2] + std::conj(t2)*prev_tmat[3];
					tmat[3] = t2*prev_tmat[2] + std::conj(t1)*prev_tmat[3];

					prev_tmat = tmat;
				}
				iso_filter->at(ix+Nx*iy) =
					(tmat[0]-tmat[1]*tmat[2]/tmat[3]) *
					std::exp(std::complex<double>(0,kN*z_foc));
			}
			else 
				iso_filter->at(ix+Nx*iy) = 0;
		}
	}
	return iso_filter;
}
