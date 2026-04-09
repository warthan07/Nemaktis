#include "Illuminator.h"

Illuminator::Illuminator(double* wavelengths, int n_wavelengths) :
    wavelengths(wavelengths, wavelengths + n_wavelengths) {}
