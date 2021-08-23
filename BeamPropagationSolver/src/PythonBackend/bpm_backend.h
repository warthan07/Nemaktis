#include <string>
#include <complex>

void run_backend_without_mask(
		std::string json_str,
		double* lc_field_vals, int n_lc_vals,
		std::complex<double>* E_field_vals, int n_E_vals);

void run_backend_with_mask(
		std::string json_str, std::string mask_formula,
		double* lc_field_vals, int n_lc_vals,
		double* mask_vals, int n_mask_vals,
		std::complex<double>* E_field_vals, int n_E_vals);
