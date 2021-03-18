#ifndef BART_SRC_QUADRATURE_QUADRATURE_SET_I_HPP_
#define BART_SRC_QUADRATURE_QUADRATURE_SET_I_HPP_

#include "quadrature/quadrature_point_i.hpp"
#include "quadrature/quadrature_types.h"
#include "problem/parameter_types.hpp"

#include <memory>
#include <optional>
#include <set>

namespace bart::quadrature {

/*! \brief Interface for quadrature sets.
 *
 * A quadrature set is a wrapper for a mapping of quadrature points
 * to an index that uniquely defines each of them. Most memeber functions are
 * for accessing points by index or vice versa. The base interface allows for
 * customization of the types of quadrature points.
 *
 * @tparam dim spatial dimension
 */
template <int dim>
class QuadratureSetI {
 public:
  using Iterator = typename std::set<std::shared_ptr<QuadraturePointI<dim>>>::iterator;
  using ConstIterator = typename std::set<std::shared_ptr<QuadraturePointI<dim>>>::const_iterator;
  using QuadraturePoint = QuadraturePointI<dim>;
  virtual ~QuadratureSetI() = default;

  /*! \brief Add a point to the quadrature set.
      * @return bool indicating if an insertion was made.
   */
  virtual auto AddPoint(std::shared_ptr<QuadraturePoint>) -> bool = 0;

  /*! \brief Returns the quadrature point from a reflection against a boundary
   * This only supports square meshes.
   */
  virtual auto GetBoundaryReflection(const std::shared_ptr<QuadraturePoint>&,
                                     problem::Boundary) const -> std::shared_ptr<QuadraturePoint> = 0;

  /// \brief Sets two points as reflections of each other.
  virtual auto SetReflection(std::shared_ptr<QuadraturePoint>, std::shared_ptr<QuadraturePoint>) -> void = 0;

  /*! \brief Returns the reflection for the provided quadrature point.
   * @return pointer to the reflection.
   */
  virtual auto GetReflection(std::shared_ptr<QuadraturePoint>) const -> std::shared_ptr<QuadraturePoint> = 0;

  /*! \brief Return the index of the provided quadrature point.
   * @return index of the refection or an empty optional if none exists.
   */
  virtual auto GetReflectionIndex(std::shared_ptr<QuadraturePoint>) const -> std::optional<int> = 0;

  /*! \brief Get quadrature point based on index.
   *
   * @param index index of quadrature point to retrieve.
   * @return pointer to quadrature point.
   */
  virtual auto GetQuadraturePoint(QuadraturePointIndex) const -> std::shared_ptr<QuadraturePoint> = 0;
  /*! \brief Get quadrature point index based on quadrature point.
   *
   * @param quadrature_point the quadrature point
   * @return int index of the point provided.
   */
  virtual auto GetQuadraturePointIndex(std::shared_ptr<QuadraturePoint>) const -> int= 0;

  /*! \brief Return the indices of the quadrature points.
   *
   * @return set containing the quadrature point indices.
   */
  virtual auto quadrature_point_indices() const -> std::set<int> = 0;

  /// \brief Returns iterator to beginning of the set of quadrature points.
  virtual auto begin() -> Iterator = 0;
  /// \brief Returns past-the-end iterator for the set of quadrature points.
  virtual auto end() -> Iterator = 0;
  /// \brief Returns const iterator to beginning of the set of quadrature points.
  virtual auto cbegin() const -> ConstIterator = 0;
  /// \brief Returns const past-the-end iterator for the set of quadrature points.
  virtual auto cend() const -> ConstIterator = 0;
  /// \brief Returns the size of the quadrature set.
  virtual auto size() const -> std::size_t = 0;
};

} // namespace bart::quadrature

#endif //BART_SRC_QUADRATURE_QUADRATURE_SET_I_HPP_
