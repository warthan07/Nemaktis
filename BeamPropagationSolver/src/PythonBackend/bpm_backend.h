#include <string>
#include <complex>

void run_backend_without_mask(
		std::string json_str,
		double* director_vals, int n_director,
		std::complex<double>* fields_vals, int n_fields);

void run_backend_with_mask(
		std::string json_str,
		double* director_vals, int n_director,
		double* mask_vals, int n_mask,
		std::complex<double>* fields_vals, int n_fields);
