#include <cstring>

#include <boost/filesystem.hpp>

#include <vtkDataArraySelection.h>
#include <vtkPointData.h>

#include "VTIReader.h"
#include "error.h"

VTIReader::VTIReader(
			std::string &filename, unsigned int lc_dim) :
		lc_dim(lc_dim) {

	Assert(
		lc_dim==3 || lc_dim==6, "Wrong field dimension for the LC solution");

	// We check if the given file exists
	boost::filesystem::path path(filename);
	if(!boost::filesystem::exists(path)) 
		throw std::string(
			"The given input vti file does not exists");

	vtkObject::GlobalWarningDisplayOff();

	reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
	reader->SetFileName(filename.c_str());
	reader->UpdateInformation();

	found_n_data = false;
	found_q_data = false;
	found_S_data = false;
	for(int i=0; i<reader->GetNumberOfPointArrays(); i++) {
		auto array_name = reader->GetPointArrayName(i);
		if(!std::strcmp(array_name, "n")) {
			found_n_data = true;
			reader->SetPointArrayStatus(array_name, 1);
		}
		else if(!std::strcmp(array_name, "Q")) {
			found_q_data = true;
			reader->SetPointArrayStatus(array_name, 1);
		}
		else if(!std::strcmp(array_name, "S")) {
			found_S_data = true;
			reader->SetPointArrayStatus(array_name, 1);
		}
		else
			reader->SetPointArrayStatus(array_name, 0);
	}
	if(lc_dim==3 && !found_n_data)
		throw std::string(
			"Could not find the director array \"n\" in the given vti file");
	if(lc_dim==6) {
		if(found_q_data) {
			// If we are in Q-tensor mode and found raw q data in the vti
			// file, we don't need n and S.
			reader->SetPointArrayStatus("n", 0);
			reader->SetPointArrayStatus("1-S/Seq", 0);
		}
		else if(!found_n_data) // If not, we need to use the n data to reconstruct q
			throw std::string(
				"Could not find the director array \"n\" in the given vti file");
	}

	reader->Update();
	output = reader->GetOutput();

	if(found_n_data)
		n_data = vtkDoubleArray::FastDownCast(
			output->GetPointData()->GetAbstractArray("n"));
	if(found_q_data)
		q_data = vtkDoubleArray::FastDownCast(
			output->GetPointData()->GetAbstractArray("Q"));
	if(found_S_data)
		S_data = vtkDoubleArray::FastDownCast(
			output->GetPointData()->GetAbstractArray("S"));
}

void VTIReader::fill_solution_vector(
		std::shared_ptr<VectorField<double> > &lc_sol,
		BasisConvention basis_convention) const {

	std::cout << "Filling data structures..." << std::endl;

	int dim[3];
	output->GetDimensions(dim);

	double spacing[3];
	output->GetSpacing(spacing);

	unsigned int I1, I2, I3;
	unsigned int Np = 0;
	switch(basis_convention) {
	case BasisConvention::XYZ:
		I1 = 0;		I2 = 1;		I3 = 2;
		break;
	case BasisConvention::YZX:
		I1 = 1;		I2 = 2;		I3 = 0;
		double norm;
		do {
			norm = std::sqrt(
				std::pow(n_data->GetComponent(Np,0),2.)+
				std::pow(n_data->GetComponent(Np,1),2.)+
				std::pow(n_data->GetComponent(Np,2),2.));
			Np++;
		} while(norm==0);
		Np--;
		break;
	case BasisConvention::ZXY:
		I1 = 2;		I2 = 0;		I3 = 1;
	}

	CartesianMesh mesh(
		{spacing[I1], spacing[I2], spacing[I3]},
		{(unsigned int)dim[I1], (unsigned int)dim[I2], (unsigned int)dim[I3]-2*Np});
	lc_sol = std::make_shared<VectorField<double> >(mesh, lc_dim);

	if(lc_dim==3) {
		for(unsigned int ix=0; ix<dim[0]-2*Np; ix++) {
			for(unsigned int iy=0; iy<dim[1]; iy++) {
				for(unsigned int iz=0; iz<dim[2]; iz++) {
					unsigned int idx[3] = {ix, iy, iz};
					unsigned int i_out = idx[I1]+dim[I1]*(idx[I2]+dim[I2]*idx[I3]);
					unsigned int i_in = Np+idx[0]+dim[0]*(idx[1]+dim[1]*idx[2]);

					if(std::isnan(n_data->GetComponent(i_in, I1)))
						throw std::string("Error: NaN value in the given vti file");
					(*lc_sol)(i_out,0) = n_data->GetComponent(i_in, I1);

					if(std::isnan(n_data->GetComponent(i_in, I2)))
						throw std::string("Error: NaN value in the given vti file");
					(*lc_sol)(i_out,1) = n_data->GetComponent(i_in, I2);

					if(std::isnan(n_data->GetComponent(i_in, I3)))
						throw std::string("Error: NaN value in the given vti file");
					(*lc_sol)(i_out,2) = n_data->GetComponent(i_in, I3);
				}
			}
		}
	}
	else if(lc_dim==6) {
		if(found_q_data) {
			for(unsigned int ix=0; ix<dim[0]; ix++) {
				for(unsigned int iy=0; iy<dim[1]; iy++) {
					for(unsigned int iz=0; iz<dim[2]; iz++) {
						unsigned int idx[3] = {ix, iy, iz};
						unsigned int i_out = idx[I1]+dim[I1]*(idx[I2]+dim[I2]*idx[I3]);
						unsigned int i_in = idx[0]+dim[0]*(idx[1]+dim[1]*idx[2]);

						for(unsigned int comp=0; comp<6; comp++) {
							Assert(
								!std::isnan(q_data->GetComponent(i_in, comp)),
								"Error: NaN value in the given vti file");
							(*lc_sol)(i_out,comp) = q_data->GetComponent(i_in, comp);
						}
					}
				}
			}
		}
		else {
			if(found_S_data) {
				#pragma omp parallel for
				for(unsigned int ix=0; ix<dim[0]; ix++) {
					for(unsigned int iy=0; iy<dim[1]; iy++) {
						for(unsigned int iz=0; iz<dim[2]; iz++) {
							unsigned int idx[3] = {ix, iy, iz};
							unsigned int i_out = idx[I1]+dim[I1]*(idx[I2]+dim[I2]*idx[I3]);
							unsigned int i_in = idx[0]+dim[0]*(idx[1]+dim[1]*idx[2]);

							(*lc_sol)(i_out,0) =
								0.5*(1-S_data->GetComponent(i_in,0))*(
									3*n_data->GetComponent(i_in,0)*n_data->GetComponent(i_in,0)-1);
							(*lc_sol)(i_out,1) =
								0.5*(1-S_data->GetComponent(i_in,0))*(
									3*n_data->GetComponent(i_in,1)*n_data->GetComponent(i_in,1)-1);
							(*lc_sol)(i_out,2) =
								0.5*(1-S_data->GetComponent(i_in,0))*(
									3*n_data->GetComponent(i_in,2)*n_data->GetComponent(i_in,2)-1);
							(*lc_sol)(i_out,3) =
								1.5*(1-S_data->GetComponent(i_in,0))*
									n_data->GetComponent(i_in,0)*n_data->GetComponent(i_in,1);
							(*lc_sol)(i_out,4) =
								1.5*(1-S_data->GetComponent(i_in,0))*
									n_data->GetComponent(i_in,0)*n_data->GetComponent(i_in,2);
							(*lc_sol)(i_out,5) =
								1.5*(1-S_data->GetComponent(i_in,0))*
									n_data->GetComponent(i_in,1)*n_data->GetComponent(i_in,2);
						}
					}
				}
			}
			else {
				#pragma omp parallel for
				for(unsigned int ix=0; ix<dim[0]; ix++) {
					for(unsigned int iy=0; iy<dim[1]; iy++) {
						for(unsigned int iz=0; iz<dim[2]; iz++) {
							unsigned int idx[3] = {ix, iy, iz};
							unsigned int i_out = idx[I1]+dim[I1]*(idx[I2]+dim[I2]*idx[I3]);
							unsigned int i_in = idx[0]+dim[0]*(idx[1]+dim[1]*idx[2]);

							(*lc_sol)(i_out,0) =
								0.5*(3*n_data->GetComponent(i_in,0)*n_data->GetComponent(i_in,0)-1);
							(*lc_sol)(i_out,1) =
								0.5*(3*n_data->GetComponent(i_in,1)*n_data->GetComponent(i_in,1)-1);
							(*lc_sol)(i_out,2) =
								0.5*(3*n_data->GetComponent(i_in,2)*n_data->GetComponent(i_in,2)-1);
							(*lc_sol)(i_out,3) =
								1.5*n_data->GetComponent(i_in,0)*n_data->GetComponent(i_in,1);
							(*lc_sol)(i_out,4) =
								1.5*n_data->GetComponent(i_in,0)*n_data->GetComponent(i_in,2);
							(*lc_sol)(i_out,5) =
								1.5*n_data->GetComponent(i_in,1)*n_data->GetComponent(i_in,2);
						}
					}
				}
			}
		}
	}
}
