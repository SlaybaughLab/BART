#ifndef BART_SRC_QUADRATURE_UTILITY_QUADRATURE_UTILITIES_H_
#define BART_SRC_QUADRATURE_UTILITY_QUADRATURE_UTILITIES_H_

#include "quadrature/ordinate_i.h"

namespace bart {

namespace quadrature {

namespace utility {

template <int dim>
std::array<double, dim> ReflectAcrossOrigin(const OrdinateI<dim>& ordinate);

template <int dim>
std::vector<std::pair<CartesianPosition<dim>, Weight>> GenerateAllPositiveX(
    const std::vector<std::pair<CartesianPosition<dim>, Weight>>&);

} // namespace utility

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_UTILITY_QUADRATURE_UTILITIES_H_
