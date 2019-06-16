#ifndef CUBICINTERPOLATEDMAPPING_H
#define CUBICINTERPOLATEDMAPPING_H

#include "Mapping.h"

template <int dim1, int dim2, typename T>
class CubicInterpolatedMapping : public InterpolatedMapping<dim1,dim2,T> {
public:
	CubicInterpolatedMapping(
		const std::shared_ptr<std::vector<Vector<dim2,T> > > &values,
		const std::shared_ptr<CartesianMesh<dim1> > &mesh,
		const std::shared_ptr<DefinitionDomain<dim1> > &def_domain);

	virtual std::shared_ptr<Mapping<dim1,dim2,T> > clone() const;

private:
	virtual void assemble_pol_weights(MultiDimIndex<dim1> &cell_origin_indices);
};

#endif
