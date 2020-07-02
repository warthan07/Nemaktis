%module "bpm_backend"

%{
#define SWIG_FILE_WITH_INIT
#include "bpm_backend.h"
%}

%include <std_string.i>
%include "numpy.i"

%init %{
    import_array();
%}

%apply (double* IN_ARRAY1, int DIM1){(double* director_vals, int n_director)};
%apply (double* IN_ARRAY1, int DIM1){(double* mask_vals, int n_mask)};
%apply (std::complex<double>* ARGOUT_ARRAY1, int DIM1){(std::complex<double>* fields_vals, int n_fields)};

%include "bpm_backend.h"
