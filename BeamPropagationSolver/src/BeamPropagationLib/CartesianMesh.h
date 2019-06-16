#ifndef CARTESIANMESH_H
#define CARTESIANMESH_H

struct CartesianMesh {
	CartesianMesh(const double (&delta)[3], const unsigned int (&N)[3]) :
		delta_x(delta[0]),
		delta_y(delta[1]),
		delta_z(delta[2]),
		Nx(N[0]),
		Ny(N[1]),
		Nz(N[2]) {}

	/**
	 * Mesh spacings.
	 */
	double delta_x, delta_y, delta_z;

	/**
	 * Number of points in each space direction.
	 */
	unsigned int Nx, Ny, Nz;
};

#endif
