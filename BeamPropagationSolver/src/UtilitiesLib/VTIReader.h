#ifndef VTIREADER_H
#define VTIREADER_H

#include <string>
#include <fstream>
#include <vector>

#include <vtkXMLImageDataReader.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>

#include "VectorField.h"
#include "PhysicsSettings.h"

class VTIReader {
public:
	VTIReader(
		std::string &filename, int lc_dim);

	void fill_solution_vector(
		std::shared_ptr<VectorField<double> > &lc_sol,
		BasisConvention basis_convention) const;

private:
	const int lc_dim;

	vtkSmartPointer<vtkXMLImageDataReader> reader;
	vtkSmartPointer<vtkImageData> output;

	vtkSmartPointer<vtkDoubleArray> n_data;
	vtkSmartPointer<vtkDoubleArray> S_data;
	vtkSmartPointer<vtkDoubleArray> q_data;

	bool found_n_data, found_q_data, found_S_data;
};

#endif
