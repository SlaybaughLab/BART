#ifndef BART_SRC_QUADRATURE_UTILITY_QUADRATURE_UTILITIES_H_
#define BART_SRC_QUADRATURE_UTILITY_QUADRATURE_UTILITIES_H_

#include "quadrature/ordinate_i.h"

namespace bart {

namespace quadrature {

namespace utility {

template <int dim>
std::array<double, dim> Reflect(const OrdinateI<dim>& ordinate);

} // namespace utility

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_UTILITY_QUADRATURE_UTILITIES_H_
