#include "quadrature/utility/quadrature_utilities.h"

#include <functional>

namespace bart {

namespace quadrature {

namespace utility {

template <int dim>
std::array<double, dim> ReflectAcrossOrigin(const OrdinateI<dim>& ordinate) {

  auto position = ordinate.cartesian_position();

  std::transform(position.begin(), position.end(),
                 position.begin(), [](double val){return -val;});

  return position;
}


template std::array<double, 1> ReflectAcrossOrigin<1>(const OrdinateI<1>&);
template std::array<double, 2> ReflectAcrossOrigin<2>(const OrdinateI<2>&);
template std::array<double, 3> ReflectAcrossOrigin<3>(const OrdinateI<3>&);

} // namespace utility

} // namespace quadrature

} // namespace bart