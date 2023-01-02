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
	coefs(coefs),
	iso_layer_thickness(settings.postprocessor.micrograph_output.iso_layer_thickness),
	iso_layer_index(settings.postprocessor.micrograph_output.iso_layer_index),
	z_foc(settings.postprocessor.micrograph_output.focalisation_z_shift),
	numerical_aperture(settings.postprocessor.micrograph_output.numerical_aperture),
	wavelengths(coefs.wavelengths()),
	q_vals(coefs.q_vals()) {

	if(iso_layer_thickness.size()!=iso_layer_index.size())
		throw std::string(
			"Inconsistent number of isotropic layers");

	// new entry for the refractive index of the medium between the last iso layer and the
	// objective, will be overwritten later on for each wavelength
	double n_out = coefs.get_nout(0.6);
	iso_layer_index.push_back(n_out);

	for(int i=0; i<iso_layer_thickness.size(); i++)
		z_foc -= iso_layer_thickness[i] * n_out / iso_layer_index[i];
}

void ScreenOutput::apply(ScreenOpticalFieldCollection &screen_optical_fields) {

	std::cout << "Calculating output transverse fields..." << std::endl;

	auto mesh = screen_optical_fields.mesh;
	Nx = mesh.Nx;
	Ny = mesh.Ny;
	delta_x = mesh.delta_x;
	delta_y = mesh.delta_y;

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
	std::vector<vtkSmartPointer<vtkDoubleArray> > vti_data_arrays(
		4*q_vals.size()*wavelengths.size());

	std::cout <<
		"    Filtering optical fields..." << std::endl;
	std::string pol_strs[2] ={"inputX_", "inputY_"};

	#pragma omp parallel for
	for(int wave_idx=0; wave_idx<wavelengths.size(); wave_idx++) {
		// FFTW arrays for this thread
		fftw_complex* field[2];
		fftw_complex* fft_field[2];
		for(int comp=0; comp<2; comp++) {
			field[comp] = static_cast<fftw_complex*>(
				fftw_malloc(sizeof(fftw_complex) * Nx*Ny));
			fft_field[comp] = static_cast<fftw_complex*>(
				fftw_malloc(sizeof(fftw_complex) * Nx*Ny));
		}

		// We override the last value of refractive index based on n_out, since it can
		// depends on the wavelength
		double wavelength = wavelengths[wave_idx];
		iso_layer_index.back() = coefs.get_nout(wavelength);

		for(int q_idx=0; q_idx<q_vals.size(); q_idx++) {
			// Propagation filter for this wavelength and incoming wavevector
			auto fourier_filter = assemble_fourier_filter(wavelength, q_vals[q_idx]);
	
			for(int pol_idx=0; pol_idx<2; pol_idx++) {
				auto real_data = vtkSmartPointer<vtkDoubleArray>::New();
				real_data->SetNumberOfComponents(3);
				real_data->SetNumberOfTuples(Nx*Ny);
				auto varname =
					"E_real_"+pol_strs[pol_idx]+
					std::to_string(wave_idx)+"_"+std::to_string(q_idx);
				real_data->SetName(varname.c_str());
				
				auto imag_data = vtkSmartPointer<vtkDoubleArray>::New();
				imag_data->SetNumberOfComponents(3);
				imag_data->SetNumberOfTuples(Nx*Ny);
				varname =
					"E_imag_"+pol_strs[pol_idx]+
					std::to_string(wave_idx)+"_"+std::to_string(q_idx);
				imag_data->SetName(varname.c_str());
	
				vti_data_arrays[4*(q_idx+q_vals.size()*wave_idx)+2*pol_idx] = real_data;
				vti_data_arrays[4*(q_idx+q_vals.size()*wave_idx)+2*pol_idx+1] = imag_data;
	
				// First, we fill the fftw input arrays with our data and switch to
				// Fourier space
				for(int comp=0; comp<2; comp++) {
					for(int i=0; i<Nx*Ny; i++) {
						field[comp][i][0] = 
							std::real(screen_optical_fields(wave_idx,q_idx,pol_idx)(i,comp));
						field[comp][i][1] = 
							std::imag(screen_optical_fields(wave_idx,q_idx,pol_idx)(i,comp));
					}
					fftw_execute_dft(forward_plan, field[comp], fft_field[comp]);
				}

				// We apply the fourier filter
				for(int i=0; i<Nx*Ny; i++) {
					auto tmp =
						fourier_filter->at(i)*
						Eigen::Vector2cd(
							std::complex<double>(fft_field[0][i][0], fft_field[0][i][1]),
							std::complex<double>(fft_field[1][i][0], fft_field[1][i][1]));
					for(int comp=0; comp<2; comp++) {
						fft_field[comp][i][0] = std::real(tmp[comp]);
						fft_field[comp][i][1] = std::imag(tmp[comp]);
					}
				}

				// We come back to real space and save the transformed data
				for(int comp=0; comp<2; comp++) {
					fftw_execute_dft(backward_plan, fft_field[comp], field[comp]);
					for(int i=0; i<Nx*Ny; i++) {
						real_data->SetComponent(i, comp, field[comp][i][0]/(Nx*Ny));
						imag_data->SetComponent(i, comp, field[comp][i][1]/(Nx*Ny));
					}
				}
				for(int i=0; i<Nx*Ny; i++) {
					real_data->SetComponent(i, 2, 0);
					imag_data->SetComponent(i, 2, 0);
				}
			}
		}

		// We free the FFTW pointers
		for(int comp=0; comp<2; comp++) {
			fftw_free(field[comp]);
			fftw_free(fft_field[comp]);
		}
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


	auto wavelength_data = vtkSmartPointer<vtkDoubleArray>::New();
	wavelength_data->SetNumberOfComponents(1);
	wavelength_data->SetNumberOfTuples(wavelengths.size());
	wavelength_data->SetName("lambda");
	for(int wave_idx=0; wave_idx<wavelengths.size(); wave_idx++)
		wavelength_data->SetComponent(wave_idx, 0, wavelengths[wave_idx]);
	vti_data->GetFieldData()->AddArray(wavelength_data);

	auto qx_data = vtkSmartPointer<vtkDoubleArray>::New();
	qx_data->SetNumberOfComponents(1);
	qx_data->SetNumberOfTuples(q_vals.size());
	qx_data->SetName("qx");
	for(unsigned int q_idx=0; q_idx<q_vals.size(); q_idx++)
		qx_data->SetComponent(q_idx, 0, q_vals[q_idx].first);
	vti_data->GetFieldData()->AddArray(qx_data);

	auto qy_data = vtkSmartPointer<vtkDoubleArray>::New();
	qy_data->SetNumberOfComponents(1);
	qy_data->SetNumberOfTuples(q_vals.size());
	qy_data->SetName("qy");
	for(unsigned int q_idx=0; q_idx<q_vals.size(); q_idx++)
		qy_data->SetComponent(q_idx, 0, q_vals[q_idx].second);
	vti_data->GetFieldData()->AddArray(qy_data);

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
		ScreenOpticalFieldCollection &screen_optical_fields,
		std::complex<double>* fields_vals) {

	std::cout << "Calculating output transverse fields..." << std::endl;

	auto mesh = screen_optical_fields.mesh;
	Nx = mesh.Nx;
	Ny = mesh.Ny;
	delta_x = mesh.delta_x;
	delta_y = mesh.delta_y;

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
		fftw_complex* field[2];
		fftw_complex* fft_field[2];
		for(int comp=0; comp<2; comp++) {
			field[comp] = static_cast<fftw_complex*>(
				fftw_malloc(sizeof(fftw_complex) * Nx*Ny));
			fft_field[comp] = static_cast<fftw_complex*>(
				fftw_malloc(sizeof(fftw_complex) * Nx*Ny));
		}

		// We override the last value of refractive index based on n_out, since it can
		// depends on the wavelength
		double wavelength = wavelengths[wave_idx];
		iso_layer_index.back() = coefs.get_nout(wavelength);

		for(int q_idx=0; q_idx<q_vals.size(); q_idx++) {
			// Propagation filter for this wavelength and incoming wavevector
			auto fourier_filter = assemble_fourier_filter(wavelengths[wave_idx], q_vals[q_idx]);

			for(int pol_idx=0; pol_idx<2; pol_idx++) {
				// First, we fill the fftw input arrays with our data and switch to
				// Fourier space
				for(int comp=0; comp<2; comp++) {
					for(int i=0; i<Nx*Ny; i++) {
						field[comp][i][0] = 
							std::real(screen_optical_fields(wave_idx,q_idx,pol_idx)(i,comp));
						field[comp][i][1] = 
							std::imag(screen_optical_fields(wave_idx,q_idx,pol_idx)(i,comp));
					}
					fftw_execute_dft(forward_plan, field[comp], fft_field[comp]);
				}

				// We apply the fourier filter
				for(int i=0; i<Nx*Ny; i++) {
					Eigen::Vector2cd tmp =
						fourier_filter->at(i)*
						Eigen::Vector2cd(
							std::complex<double>(fft_field[0][i][0], fft_field[0][i][1]),
							std::complex<double>(fft_field[1][i][0], fft_field[1][i][1]));
					for(int comp=0; comp<2; comp++) {
						fft_field[comp][i][0] = std::real(tmp[comp]);
						fft_field[comp][i][1] = std::imag(tmp[comp]);
					}
				}

				// We come back to real space and save the transformed data
				for(int comp=0; comp<2; comp++) {
					fftw_execute_dft(backward_plan, fft_field[comp], field[comp]);
					for(int i=0; i<Nx*Ny; i++)
						fields_vals[i+Nx*Ny*(comp+2*pol_idx+4*(q_idx+q_vals.size()*wave_idx))]=
							std::complex<double>(field[comp][i][0],field[comp][i][1]) /
							double(Nx*Ny);
				}
			}
		}

		// We free the FFTW pointers
		for(int comp=0; comp<2; comp++) {
			fftw_free(field[comp]);
			fftw_free(fft_field[comp]);
		}
	}

	// We free the FFTW pointers
	fftw_destroy_plan(forward_plan);
	fftw_destroy_plan(backward_plan);
	fftw_free(in);
	fftw_free(out);
}

std::shared_ptr<std::vector<Eigen::Matrix2cd> > ScreenOutput::assemble_fourier_filter(
		double wavelength, std::pair<double,double> q_val) const {

	auto fourier_filter = std::make_shared<std::vector<Eigen::Matrix2cd> >(Nx*Ny);
	double k0 = 2*PI/wavelength;
	double n_out = iso_layer_index.back();

	// total and local transfer matrix for forward and backward modes 
	Eigen::Matrix4cd tmat, layer_tmat; 
	// temporary matrices for the assembly of the transfer matrix
	Eigen::Matrix2cd tmp; 
	Eigen::Matrix2d index_mat, next_inv_index_mat, identity, tp, tm;
	identity.setIdentity();

	#pragma omp parallel for firstprivate(tmat, layer_tmat, \
		index_mat, next_inv_index_mat, identity)
	for(int iy=0; iy<Ny; iy++) {
		for(int ix=0; ix<Nx; ix++) {
			double kx = (ix<double(Nx/2.)) ?
				2*PI*ix/(delta_x*Nx) :
				-2*PI*(Nx-ix)/(delta_x*Nx);
			double ky = (iy<double(Ny/2.)) ?
				2*PI*iy/(delta_y*Ny) :
				-2*PI*(Ny-iy)/(delta_y*Ny);
			Eigen::Vector2d q_vec = {kx/k0+q_val.first, ky/k0+q_val.second};
			double q_norm_sqr = std::pow(q_vec.norm(),2);
			double q0_sqr = std::sqrt(std::pow(q_val.first,2)+std::pow(q_val.second,2));

			if(q_vec.norm()<numerical_aperture) {
				tmat.setIdentity();

				for(int p=0; p<iso_layer_thickness.size(); p++) {
					double np = iso_layer_index[p];
					double next_np = iso_layer_index[p+1];
					double neff = std::sqrt(np*np-q_norm_sqr);
					double next_neff = std::sqrt(next_np*next_np-q_norm_sqr);

					index_mat <<
						neff+q_vec[0]*q_vec[0]/neff,
						q_vec[0]*q_vec[1]/neff,
						q_vec[0]*q_vec[1]/neff,
						neff+q_vec[1]*q_vec[1]/neff;
					next_inv_index_mat <<
						(1-std::pow(q_vec[0]/next_np,2))/next_neff,
						-q_vec[0]*q_vec[1]/(next_np*next_np*next_neff),
						-q_vec[0]*q_vec[1]/(next_np*next_np*next_neff),
						(1-std::pow(q_vec[1]/next_np,2))/next_neff;

					tp = 0.5*(identity+next_inv_index_mat*index_mat);
					tm = identity-tp;

					auto g = std::exp(std::complex<double>(0,
						(neff-std::sqrt(np*np-q0_sqr))*k0*iso_layer_thickness[p]));
					layer_tmat.topLeftCorner<2,2>() = tp*g;
					layer_tmat.topRightCorner<2,2>() = tm*std::conj(g);
					layer_tmat.bottomLeftCorner<2,2>() = tm*g;
					layer_tmat.bottomRightCorner<2,2>() = tp*std::conj(g);

					tmat.noalias() = layer_tmat*tmat;
				}

				// We invert the (1,1) block and calculate the reduced transfer matrix for
				// forward modes
				auto det = tmat(2,2)*tmat(3,3)-tmat(2,3)*tmat(3,2);
				tmp << tmat(3,3)/det, -tmat(2,3)/det, -tmat(3,2)/det, tmat(2,2)/det;
				fourier_filter->at(ix+Nx*iy) =
					tmat.topLeftCorner<2,2>() -
					tmat.topRightCorner<2,2>()*tmp*tmat.bottomLeftCorner<2,2>();
				fourier_filter->at(ix+Nx*iy) *= std::exp(std::complex<double>(0,
					(std::sqrt(n_out*n_out-q_norm_sqr)-std::sqrt(n_out*n_out-q0_sqr))*k0*z_foc));
			}
			else 
				fourier_filter->at(ix+Nx*iy).setZero();
		}
	}
	return fourier_filter;
}
