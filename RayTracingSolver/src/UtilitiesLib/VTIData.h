#ifndef VTIDATA_H
#define VTIDATA_H

#include <string>
#include <vector>

#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkSmartPointer.h>


class RayData {
public:
	RayData(
		unsigned long Nx, unsigned long Ny, unsigned long Nz,
		double dx, double dy, double dz);

	void write(std::string basedir, std::string filename);

	vtkSmartPointer<vtkImageData> vti_data;
	vtkSmartPointer<vtkDoubleArray> deflection, opt_length, moment; 
	vtkSmartPointer<vtkDoubleArray> ampl[2];
};


class FieldsData {
public:
	FieldsData(
		unsigned long Nx, unsigned long Ny, unsigned long Nz,
		double dx, double dy, double dz, 
		std::vector<double> &wavelengths);

	void write(std::string basedir, std::string filename);

	vtkSmartPointer<vtkImageData> vti_data;
	std::vector<vtkSmartPointer<vtkDoubleArray> > E_real[2], E_imag[2];
	vtkSmartPointer<vtkDoubleArray> ray_multiplicity;
};


#endif
