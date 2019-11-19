#ifndef BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_TYPES_H_
#define BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_TYPES_H_

#include <array>
#include <utility>

namespace bart {

namespace quadrature {

namespace angular {

/*! Weights for ordinates */
using Weight [[deprecated("Replaced by named type Weight")]] = double;

/*! Ordinate points for quadrature sets */
template <int dim>
using Ordinate [[deprecated("Replaced by class OrdinateI")]] = std::array<double, dim>;

/*! Quadrature points */
template <int dim>
using QuadraturePoint [[deprecated("Replaced by class QuadraturePointI")]] = std::pair<Weight, Ordinate<dim>>;

} // namespace angular

} // namespace quadrature

} // namespace bart

#endif // BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_TYPES_H_