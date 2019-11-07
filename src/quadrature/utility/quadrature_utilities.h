#ifndef BART_SRC_QUADRATURE_UTILITY_QUADRATURE_UTILITIES_H_
#define BART_SRC_QUADRATURE_UTILITY_QUADRATURE_UTILITIES_H_

#include "quadrature/quadrature_point_i.h"
#include "quadrature/ordinate_i.h"

namespace bart {

namespace quadrature {

namespace utility {

template <int dim>
std::array<double, dim> ReflectAcrossOrigin(const OrdinateI<dim>& ordinate);

template <int dim>
std::vector<std::pair<CartesianPosition<dim>, Weight>> GenerateAllPositiveX(
    const std::vector<std::pair<CartesianPosition<dim>, Weight>>&);

template <int dim>
struct quadrature_point_compare{
  bool operator() (const std::shared_ptr<quadrature::QuadraturePointI<dim>>& lhs,
                   const std::shared_ptr<quadrature::QuadraturePointI<dim>>& rhs) {
    return lhs->cartesian_position() < rhs->cartesian_position();
  }
};

} // namespace utility

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_UTILITY_QUADRATURE_UTILITIES_H_
