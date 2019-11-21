#ifndef BART_SRC_QUADRATURE_QUADRATURE_SET_I_H_
#define BART_SRC_QUADRATURE_QUADRATURE_SET_I_H_

#include "quadrature/quadrature_point_i.h"

#include <memory>
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

  /// \brief Sets two points as reflections of each other.
  virtual void SetReflection(std::shared_ptr<QuadraturePointI<dim>>,
                             std::shared_ptr<QuadraturePointI<dim>>) = 0;

  /*! \brief Returns the reflection for the provided quadrature point.
   * @return pointer to the reflection.
   */
  virtual std::shared_ptr<QuadraturePointI<dim>> GetReflection(
      std::shared_ptr<QuadraturePointI<dim>>) const = 0;

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
