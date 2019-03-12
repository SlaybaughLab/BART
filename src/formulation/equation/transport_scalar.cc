#include "formulation/equation/transport_scalar.h"

namespace bart {

namespace formulation {

namespace equation {

template <int dim>
void TransportScalar<dim>::FillCellVariableBilinear(Matrix&,
                                                    const CellPtr& ,
                                                    const GroupNumber) const {}

template <int dim>
void TransportScalar<dim>::FillBoundaryVariableBilinear(
    Matrix &,
    const CellPtr &,
    const GroupNumber,
    const FaceNumber,
    const BoundaryType) const {}

template class TransportScalar<1>;
template class TransportScalar<2>;
template class TransportScalar<3>;

} // namespace equation

} // namespace formulation

} // namespace bart