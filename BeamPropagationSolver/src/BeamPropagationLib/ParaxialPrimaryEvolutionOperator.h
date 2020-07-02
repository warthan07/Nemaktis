#ifndef PARAXIALPRIMARYEVOLUTIONOPERATOR_H
#define PARAXIALPRIMARYEVOLUTIONOPERATOR_H

#include "BaseBPMOperator.h"
#include "BlockSparseMatrix.h"

class ParaxialPrimaryEvolutionOperator : public BaseBPMOperator {
public:
	ParaxialPrimaryEvolutionOperator(
		const PermittivityTensorField& eps,
		double wavelength, const RootSettings &settings);

	/**
	 * Apply this evolution operator to a transverse optical field.
	 */
	virtual int apply(TransverseOpticalField &src);
	/**
	 * Reinit this operator for the current transverse plane.
	 */
	virtual void update();

private:
	/**
	 * Reinit D for the current transverse plane.
	 */
	void update_D();
	/**
	 * Reinit D1 for the current transverse plane.
	 */
	void update_D1();

	/**
	 * Matrix of the operator \partial_x in a row-major format.
	 */
	Eigen::SparseMatrix<double,Eigen::RowMajor> d_ovr_dX_csr;
	/**
	 * Matrix of the operator \partial_x in a col-major format.
	 */
	Eigen::SparseMatrix<double,Eigen::ColMajor> d_ovr_dX_csc;
	/**
	 * Matrix of the operator \partial_y in a row-major format.
	 */
	Eigen::SparseMatrix<double,Eigen::RowMajor> d_ovr_dY_csr;
	/**
	 * Matrix of the operator \partial_y in a col-major format.
	 */
	Eigen::SparseMatrix<double,Eigen::ColMajor> d_ovr_dY_csc;

	/**
	 * Isotropic part of the diffraction matrix D, defined as the action
	 * of the differential operator
	 * [(\partial^2_x+\partial^2_y)\delta_ab - \partial_a\partial_b],
	 * with \delta the kronecker symbol and a,b=x or y.
	 */
	SplittedBlockSparseMatrix<double,Eigen::RowMajor> D0_matrix;
	/**
	 * The anisotropic part of the diffraction matrix D can be written as
	 * D1*K where K is defined below. This member attritbute stores the
	 * matrix D1, defined as the action of the differential operator
	 * [\partial_a eps_zz^-1 \partial_b] 
	 */
	SplittedBlockSparseMatrix<double,Eigen::ColMajor> D1_matrix;
	/*
	 * Diffraction matrix D = D0 + D1*K. D0 and D1 are defined above in the
	 * attributes D0_matrix and D1_matrix. The matrix K corresponds to the pointwise
	 * operator (eps_ab - exp_az*eps_zb/eps_zz), whose values can be found directly
	 * in the PermittivityTensorField object (method xx_tr, yy_tr and xy_tr). Note
	 * that this matrix and D0 are stored in RowMajor format, since it is a more
	 * efficient storage scheme for matrix-vector product, while the matrix D1 is
	 * stored in the default ColMajor format because it greatly simplifies the
	 * calculation of the product D1*K.
	 */
	SplittedBlockSparseMatrix<double,Eigen::RowMajor> D_matrix;
	/**
	 * Main matrix of this evolution operator, defined as R = i*W + (D*sqrt(K)^-1 +
	 * sqrt(K)^-1*D)/4, where D and K are defined above and W is the walkoff operator 
	 * (v_a (\partial_b .) + \partial_a(v_b .)) / 2, with v_a=eps_az/eps_zz and
	 * a,b=x,y.
	 */
	SplittedBlockSparseMatrix<std::complex<double>,Eigen::RowMajor> R_matrix;

	/**
	 * Total number of points in a transverse plane
	 */
	const int N;
	/**
	 * Number of woodbury iterations when applying the evolution operator 
	 */
	const int N_woodbury_steps;
	/**
	 * A temporary transverse field used in the apply method.
	 */
	TransverseOpticalField tmp;
	/**
	 * The evolution operator is rigourously defined as exp(2*mu*R_matrix), where
	 * mu=I*delta_Z/2 is stored in this variable.
	 */
	std::complex<double> mu;
};

#endif
