#include <boost/filesystem.hpp>
#include <fftw3.h>

#include "Simulation.h"

template <int dim>
Simulation<dim>::Simulation(
		json settings, double lc_thickness, unsigned int N_lc_steps,
		std::shared_ptr<CubicInterpolatedMapping<dim,3,double> > &n_field) :
	mat_properties(settings.at("Material properties")),
	N_lc_steps(N_lc_steps),
	lc_thickness(lc_thickness),
	n_field(n_field) {

	// We initialize the wavelength array
	auto &j = settings.at("Light source");
	if(j.find("Wavelengths") != j.end()) {
		wavelengths = j.at("Wavelengths").get<std::vector<double> >();
		N_wavelengths = wavelengths.size();
	}
	else {
		N_wavelengths = j.at("N wavelengths");
		double mean_lambda = j.at("Mean wavelength");
		double fwhm = j.at("Spectral FWHM");

		if(N_wavelengths==1)
			wavelengths.push_back(j.at("Mean wavelength"));
		else {
			double lambda;
			for(int i=0; i<N_wavelengths; i++) {
				wavelengths.push_back(mean_lambda - 0.5*fwhm + fwhm*i/(N_wavelengths-1));
			}
		}
	}

	// We initialize the ODE for the isotropic layers
	lower_iso_thickness = std::accumulate(
		mat_properties.e_lower_layers.begin(),
		mat_properties.e_lower_layers.end(), 0.);
	upper_iso_thickness = std::accumulate(
		mat_properties.e_upper_layers.begin(),
		mat_properties.e_upper_layers.end(), 0.);

	lower_iso_odes.push_back(std::make_shared<IsotropicODEFunction<dim> >(
		std::pow(mat_properties.n_out, 2.),
		std::make_shared<SlabDomain<dim> >(
			BasisVector<dim>(dim-1),
			-lower_iso_thickness-lc_thickness*1.5,
			-lower_iso_thickness-lc_thickness/2.)));

	double z_interface = -lower_iso_thickness-lc_thickness/2.;
	for(int i=0; i<mat_properties.n_lower_layers.size(); i++) {
		lower_iso_odes.push_back(std::make_shared<IsotropicODEFunction<dim> >(
			std::pow(mat_properties.n_lower_layers[i],2.),
			std::make_shared<SlabDomain<dim> >(
				BasisVector<dim>(dim-1), z_interface,
				z_interface + mat_properties.e_lower_layers[i])));
		z_interface += mat_properties.e_lower_layers[i];
	}

	z_interface = lc_thickness/2;
	for(int i=0; i<mat_properties.e_upper_layers.size(); i++) {
		upper_iso_odes.push_back(std::make_shared<IsotropicODEFunction<dim> >(
			std::pow(mat_properties.n_upper_layers[i],2.),
			std::make_shared<SlabDomain<dim> >(
				BasisVector<dim>(dim-1), z_interface,
				z_interface + mat_properties.e_upper_layers[i])));
		z_interface += mat_properties.e_upper_layers[i];
	}
	upper_iso_odes.push_back(std::make_shared<IsotropicODEFunction<dim> >(
		std::pow(mat_properties.n_out,2.),
		std::make_shared<SlabDomain<dim> >(
			BasisVector<dim>(dim-1), z_interface,
			z_interface+lc_thickness)));

	// We initalize the other parameters
	N_refinement_cycles =
		settings.at("InverseScreenMap parameters").at("N refinement cycles");
	hc_parameters = settings.at("InverseScreenMap parameters");

	viz_settings = settings.at("Visualisation");
	numerical_aperture = viz_settings.at("Screen output").at("Numerical aperture");

	bulk_basename = viz_settings.at("Bulk output").at("Base name");
	screen_basename = viz_settings.at("Screen output").at("Base name");
	basedir = viz_settings.at("Results folder name");

	bulk_ray_output = viz_settings.at("Bulk output").at("Export ray data");
	bulk_fields_output = viz_settings.at("Bulk output").at("Export reconstructed fields");
	screen_ray_output = viz_settings.at("Screen output").at("Export ray data");
	screen_fields_output = viz_settings.at("Screen output").at("Export reconstructed fields");

	auto widths = parse_Vector<dim-1,double>(viz_settings, "Target output widths");
	auto dims = parse_Vector<dim-1,unsigned long>(viz_settings, "Target N pixels per dim");
	auto origin = -widths/2;

	Vector<dim-1,unsigned long> coarse_dims;
	for(int i=0; i<dim-1; i++) {
		coarse_dims(i) =
			std::floor(1+(dims(i)-1)/std::pow(2,N_refinement_cycles));
		dims(i) = 1+(coarse_dims(i)-1)*std::pow(2,N_refinement_cycles);
	}
	coarse_horiz_mesh = std::make_shared<CartesianMesh<dim-1> >(
		origin, widths, coarse_dims, false);
	full_horiz_mesh = std::make_shared<CartesianMesh<dim-1> >(
		origin, widths, dims, false);

	switch(dim) {
		case 2:
			typical_length = widths(0);
			break;
		case 3:
			typical_length = std::max(widths(0),widths(1));
			break;
	}

	double bulk_z_length = lc_thickness;
	if(settings.at("Geometry").at("Sample type")=="Droplet")
		bulk_z_length -= 2*settings.at("Geometry").at("Droplet sample parameters").at(
			"Distance from upper sample plate").get<double>();
	z_step = bulk_z_length / (N_lc_steps-1);

	double n_above_lc =
		mat_properties.n_upper_layers.size()>0 ?
		mat_properties.n_upper_layers[0] : mat_properties.n_out;
	z_foc_1 =
		bulk_z_length*(1-n_above_lc/mat_properties.n_lc)/2. +
		(lc_thickness-bulk_z_length)*(1-n_above_lc/mat_properties.n_host)/2.;
	z_foc_2 = 
		bulk_z_length*(1-mat_properties.n_out/mat_properties.n_lc)/2. +
		(lc_thickness-bulk_z_length)*(1-mat_properties.n_out/mat_properties.n_host)/2.;
	for(int i=0; i<mat_properties.e_upper_layers.size(); i++)
		z_foc_2 +=
			mat_properties.e_upper_layers[i] *
			(1-mat_properties.n_out/mat_properties.n_upper_layers[i]);

	switch(dim) {
		case 2:
			N_horiz_pixels = dims(0);
			if(bulk_ray_output) {
				bulk_extra_data = std::make_shared<RayData>(
					1, dims(0), N_lc_steps,
					0., widths(0)/(dims(0)-1), z_step);
				bulk_ordi_data = std::make_shared<RayData>(
					1, dims(0), N_lc_steps,
					0., widths(0)/(dims(0)-1), z_step);
			}
			if(bulk_fields_output)
				bulk_fields_data = std::make_shared<FieldsData>(
					1, dims(0), N_lc_steps,
					0., widths(0)/(dims(0)-1), z_step,
					wavelengths);
			if(screen_ray_output) {
				screen_extra_data = std::make_shared<RayData>(
					1, dims(0), 1,
					0., widths(0)/(dims(0)-1), 0.);
				screen_ordi_data = std::make_shared<RayData>(
					1, dims(0), 1,
					0., widths(0)/(dims(0)-1), 0.);
			}
			if(screen_fields_output)
				screen_fields_data = std::make_shared<FieldsData>(
					1, dims(0), 1,
					0., widths(0)/(dims(0)-1), 0.,
					wavelengths);
			break;

		case 3:
			N_horiz_pixels = dims(0)*dims(1);
			if(bulk_ray_output) {
				bulk_extra_data = std::make_shared<RayData>(
					dims(0), dims(1), N_lc_steps,
					widths(0)/(dims(0)-1), widths(1)/(dims(1)-1), z_step);
				bulk_ordi_data = std::make_shared<RayData>(
					dims(0), dims(1), N_lc_steps,
					widths(0)/(dims(0)-1), widths(1)/(dims(1)-1), z_step);
			}
			if(bulk_fields_output)
				bulk_fields_data = std::make_shared<FieldsData>(
					dims(0), dims(1), N_lc_steps,
					widths(0)/(dims(0)-1), widths(1)/(dims(1)-1), z_step,
					wavelengths);
			if(screen_ray_output) {
				screen_extra_data = std::make_shared<RayData>(
					dims(0), dims(1), 1,
					widths(0)/(dims(0)-1), widths(1)/(dims(1)-1), 0.);
				screen_ordi_data = std::make_shared<RayData>(
					dims(0), dims(1), 1,
					widths(0)/(dims(0)-1), widths(1)/(dims(1)-1), 0.);
			}
			if(screen_fields_output)
				screen_fields_data = std::make_shared<FieldsData>(
					dims(0), dims(1), 1,
					widths(0)/(dims(0)-1), widths(1)/(dims(1)-1), 0.,
					wavelengths);
			break;
	}
}

template <int dim>
void Simulation<dim>::apply_fourier_iso_layer_filter() {

	int dims[3]; 		screen_fields_data->vti_data->GetDimensions(dims);
	double spacings[3];	screen_fields_data->vti_data->GetSpacing(spacings);

	fftw_complex *in = static_cast<fftw_complex*>(
		fftw_malloc(sizeof(fftw_complex) * dims[0]*dims[1]));
	fftw_complex *out = static_cast<fftw_complex*>(
		fftw_malloc(sizeof(fftw_complex) * dims[0]*dims[1]));
	fftw_plan forward_plan = (dim==3) ?
		fftw_plan_dft_2d(dims[1], dims[0], in, out, FFTW_FORWARD, FFTW_ESTIMATE) :
		fftw_plan_dft_1d(dims[1], in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan backward_plan = (dim==3) ?
		fftw_plan_dft_2d(dims[1], dims[0], in, out, FFTW_BACKWARD, FFTW_ESTIMATE) :
		fftw_plan_dft_1d(dims[1], in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

	#pragma omp parallel for
	for(unsigned int wave_idx=0; wave_idx<N_wavelengths; wave_idx++) {
		// First we assemble the filter
		auto iso_filter =
			std::make_shared<std::vector<std::complex<double> > >(dims[0]*dims[1]);
		double k0 = 2*M_PI/wavelengths[wave_idx];

		// transfer matrix in column-major order
		std::vector<std::complex<double> > tmat(4), prev_tmat(4); 
		std::complex<double> t1, t2;

		for(unsigned int iy=0; iy<dims[1]; iy++) {
			for(unsigned int ix=0; ix<dims[0]; ix++) {
				double kx = (ix<double(dims[0]/2.)) ?
					2*M_PI*ix/(spacings[0]*(dims[0]-1)) :
					-2*M_PI*(dims[0]-ix)/(spacings[0]*(dims[0]-1));
				if(dim==2)
					kx = 0;
				double ky = (iy<double(dims[1]/2.)) ?
					2*M_PI*iy/(spacings[1]*(dims[1]-1)) :
					-2*M_PI*(dims[1]-iy)/(spacings[1]*(dims[1]-1));
				double k = std::sqrt(kx*kx+ky*ky);

				if(k<k0*numerical_aperture) {
					tmat[0] = 1;
					tmat[1] = 0;
					tmat[2] = 0;
					tmat[3] = 1;

					prev_tmat = tmat;

					double kN = std::sqrt(k0*k0-k*k);

					for(unsigned int p=0; p<mat_properties.e_upper_layers.size(); p++) {
						double np = mat_properties.n_upper_layers[p];
						double ep = mat_properties.e_upper_layers[p];
						double kp = std::sqrt(std::pow(k0*np,2.)-k*k);
						double next_kp = (p+1<mat_properties.e_upper_layers.size()) ?
							std::sqrt(std::pow(k0*mat_properties.n_upper_layers[p+1],2.)-k*k) :
							std::sqrt(std::pow(k0*mat_properties.n_out,2.)-k*k);

						t1 = 0.5*(1+kp/next_kp)*std::exp(std::complex<double>(0,(kp-kN)*ep));
						t2 = 0.5*(1-kp/next_kp)*std::exp(std::complex<double>(0,(kp-kN)*ep));

						tmat[0] = t1*prev_tmat[0] + std::conj(t2)*prev_tmat[1];
						tmat[1] = t2*prev_tmat[0] + std::conj(t1)*prev_tmat[1];
						tmat[2] = t1*prev_tmat[2] + std::conj(t2)*prev_tmat[3];
						tmat[3] = t2*prev_tmat[2] + std::conj(t1)*prev_tmat[3];

						prev_tmat = tmat;
					}
					iso_filter->at(ix+dims[0]*iy) =
						(tmat[0]-tmat[1]*tmat[2]/tmat[3]) *
						std::exp(std::complex<double>(0,kN*z_foc_2));
				}
				else 
					iso_filter->at(ix+dims[0]*iy) = 0;
			}
		}

		// FFTW arrays for this thread
		fftw_complex *field = static_cast<fftw_complex*>(
			fftw_malloc(sizeof(fftw_complex) * dims[0]*dims[1]));
		fftw_complex *fft_field = static_cast<fftw_complex*>(
			fftw_malloc(sizeof(fftw_complex) * dims[0]*dims[1]));

		for(unsigned int pol_idx=0; pol_idx<2; pol_idx++) {
			for(unsigned int comp=0; comp<2; comp++) {
				// First, we fill the fftw input array with our data
				for(unsigned int i=0; i<dims[0]*dims[1]; i++) {
					field[i][0] =
						screen_fields_data->E_real[pol_idx][wave_idx]->GetComponent(i,comp);
					field[i][1] = 
						screen_fields_data->E_imag[pol_idx][wave_idx]->GetComponent(i,comp);
				}

				// We switch to fourrier space and apply the iso filter
				fftw_execute_dft(forward_plan, field, fft_field);
				for(unsigned int i=0; i<dims[0]*dims[1]; i++) {
					auto val = std::complex<double>(fft_field[i][0], fft_field[i][1]);
					fft_field[i][0] = std::real(val * iso_filter->at(i));
					fft_field[i][1] = std::imag(val * iso_filter->at(i));
				}

				// We come back to real space and save the transformed
				// data
				fftw_execute_dft(backward_plan, fft_field, field);

				for(unsigned int i=0; i<dims[0]*dims[1]; i++) {
					screen_fields_data->E_real[pol_idx][wave_idx]->SetComponent(
						i, comp, field[i][0]/(dims[0]*dims[1]));
					screen_fields_data->E_imag[pol_idx][wave_idx]->SetComponent(
						i, comp, field[i][1]/(dims[0]*dims[1]));
				}
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

template class Simulation<2>;
template class Simulation<3>;
