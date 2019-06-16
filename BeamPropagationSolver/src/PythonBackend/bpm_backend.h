#include <string>
#include <complex>

void run_backend(
		std::string json_str,
		double* director_vals, int n_director,
		std::complex<double>* fields_vals, int n_fields);
