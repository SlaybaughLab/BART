#ifndef BART_SRC_QUADRATURE_QUADRATURE_SET_I_H_
#define BART_SRC_QUADRATURE_QUADRATURE_SET_I_H_

#include "quadrature/quadrature_point_i.hpp"
#include "quadrature/quadrature_types.h"
#include "problem/parameter_types.hpp"

#include <memory>
#include <optional>
#include <set>

namespace bart {

namespace quadrature {

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
  virtual ~QuadratureSetI() = default;

  /*! \brief Add a point to the quadrature set.
      * @return bool indicating if an insertion was made.
   */
  virtual bool AddPoint(std::shared_ptr<QuadraturePointI<dim>>) = 0;

  /*! \brief Returns the quadrature point from a reflection against a boundary
   * This only supports square meshes.
   */
  virtual std::shared_ptr<QuadraturePointI<dim>> GetBoundaryReflection(
      const std::shared_ptr<QuadraturePointI<dim>>&,
      const problem::Boundary) const = 0;

  /// \brief Sets two points as reflections of each other.
  virtual void SetReflection(std::shared_ptr<QuadraturePointI<dim>>,
                             std::shared_ptr<QuadraturePointI<dim>>) = 0;

  /*! \brief Returns the reflection for the provided quadrature point.
   * @return pointer to the reflection.
   */
  virtual std::shared_ptr<QuadraturePointI<dim>> GetReflection(
      std::shared_ptr<QuadraturePointI<dim>>) const = 0;

  /*! \brief Return the index of the provided quadrature point.
   * @return index of the refection or an empty optional if none exists.
   */
  virtual std::optional<int> GetReflectionIndex(
      std::shared_ptr<QuadraturePointI<dim>>) const = 0;

  /*! \brief Get quadrature point based on index.
   *
   * @param index index of quadrature point to retrieve.
   * @return pointer to quadrature point.
   */
  virtual std::shared_ptr<QuadraturePointI<dim>> GetQuadraturePoint(
      QuadraturePointIndex index) const = 0;
  /*! \brief Get quadrature point index based on quadrature point.
   *
   * @param quadrature_point the quadrature point
   * @return int index of the point provided.
   */
  virtual int GetQuadraturePointIndex(
      std::shared_ptr<QuadraturePointI<dim>> quadrature_point) const = 0;

  /*! \brief Return the indices of the quadrature points.
   *
   * @return set containing the quadrature point indices.
   */
  virtual std::set<int> quadrature_point_indices() const = 0;

  /// \brief Returns iterator to beginning of the set of quadrature points.
  virtual Iterator begin() = 0;
  /// \brief Returns past-the-end iterator for the set of quadrature points.
  virtual Iterator end() = 0;
  /// \brief Returns const iterator to beginning of the set of quadrature points.
  virtual ConstIterator cbegin() const = 0;
  /// \brief Returns const past-the-end iterator for the set of quadrature points.
  virtual ConstIterator cend() const = 0;
  /// \brief Returns the size of the quadrature set.
  virtual std::size_t size() const = 0;
};

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_QUADRATURE_SET_I_H_
