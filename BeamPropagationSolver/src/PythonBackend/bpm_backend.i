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

%apply (double* IN_ARRAY1, int DIM1){(double* lc_field_vals, int n_lc_vals)};
%apply (double* IN_ARRAY1, int DIM1){(double* mask_vals, int n_mask_vals)};
%apply (std::complex<double>* INPLACE_ARRAY_FLAT, int DIM_FLAT){(std::complex<double>* E_field_vals, int n_E_vals)};

%include "bpm_backend.h"
