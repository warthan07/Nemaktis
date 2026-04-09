#include <vector>

class Illuminator {
public:    
    Illuminator(double* wavelengths, int n_wavelengths);

private:
    std::vector<double> wavelengths;
};