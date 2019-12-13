#ifndef CARTESIANMESH_H
#define CARTESIANMESH_H

struct CartesianMesh {
	CartesianMesh(const double (&delta)[3], const int (&N)[3]) :
		Nx(N[0]),
		Ny(N[1]),
		Nz(N[2]),
		delta_x(delta[0]),
		delta_y(delta[1]),
		delta_z(delta[2]),
		width_x(delta[0]*(N[0]-1)),
		width_y(delta[1]*(N[1]-1)),
		width_z(delta[2]*(N[2]-1)),
		origin_x(-width_x/2.),
		origin_y(-width_y/2.),
		origin_z(-width_z/2.) {}
		

	/**
	 * Number of points in each space direction.
	 */
	int Nx, Ny, Nz;
	/**
	 * Mesh spacings.
	 */
	double delta_x, delta_y, delta_z;
	/**
	 * Mesh widths.
	 */
	double width_x, width_y, width_z;
	/**
	 * Mesh origin.
	 */
	double origin_x, origin_y, origin_z;
};

#endif
