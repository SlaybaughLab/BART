#ifndef BART_SRC_QUADRATURE_UTILITY_QUADRATURE_UTILITIES_H_
#define BART_SRC_QUADRATURE_UTILITY_QUADRATURE_UTILITIES_H_

#include "quadrature/quadrature_point_i.hpp"
#include "quadrature/ordinate_i.hpp"

namespace bart {

namespace quadrature {

namespace utility {

/*! \brief Reflects a given ordinate across the origin.
 * @tparam dim spatial dimension.
 * @param ordinate ordinate to be reflected.
 * @return array containing the cartesian coordiante of the reflection.
 */
template <int dim>
std::array<double, dim> ReflectAcrossOrigin(const OrdinateI<dim>& ordinate);

/*! \brief Generates all pairs of positions and weights in the positive X quadrants.
 *
 * This function takes each pair of cartesian positions and weights and
 * generates a vector with that pair and all reflections across the Y-axis and
 * Z-axis (as appropriate for the spatial dimension). This function therefore
 * does nothing for dim = 1. This effectively takes a single-quadrant quadrature
 * set and "fills" all the positive-x quadrants. For 2D the resulting vector
 * will be twice as large and for 3D, four times as large.
  * @tparam dim spatial dimension.
 * @return a vector containing the original points and all reflections.
 */
template <int dim>
std::vector<std::pair<CartesianPosition<dim>, Weight>> GenerateAllPositiveX(
    const std::vector<std::pair<CartesianPosition<dim>, Weight>>&);

/*! \brief Struct to compare two quadrature points based on cartesian position.
 * This is required for using quadrature points in any associative container
 * such as maps, sets, etc.
 * @tparam dim spatial dimension.
 */
template <int dim>
struct quadrature_point_compare{
  /*! \brief Parenthesis operator that provides a less-than comparison.
   *
   * @param lhs left-hand-side of comparison.
   * @param rhs right-hand-side of comparison.
   * @return bool indicating if left-hand-side is less than right-hand-side
   */
  bool operator() (const std::shared_ptr<quadrature::QuadraturePointI<dim>>& lhs,
                   const std::shared_ptr<quadrature::QuadraturePointI<dim>>& rhs) {
    return lhs->cartesian_position() < rhs->cartesian_position();
  }

  /*! \brief Const parenthesis operator that provides a less-than comparison.
   *
   * @param lhs left-hand-side of comparison.
   * @param rhs right-hand-side of comparison.
   * @return bool indicating if left-hand-side is less than right-hand-side
   */
  bool operator() (const std::shared_ptr<quadrature::QuadraturePointI<dim>>& lhs,
                   const std::shared_ptr<quadrature::QuadraturePointI<dim>>& rhs) const {
    return lhs->cartesian_position() < rhs->cartesian_position();
  }
};

} // namespace utility

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_UTILITY_QUADRATURE_UTILITIES_H_
