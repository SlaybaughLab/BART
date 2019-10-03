#include "quadrature/ordinate.h"

namespace bart {

namespace quadrature {

template<int dim>
Ordinate<dim>::Ordinate(CartesianPosition<dim> position)
    : cartesian_position_(position.get())
{}

template class Ordinate<1>;
template class Ordinate<2>;
template class Ordinate<3>;

} // namespace quadrature

} // namespace bart
