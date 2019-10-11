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

  virtual bool AddPoint(std::shared_ptr<QuadraturePointI<dim>>) = 0;

  virtual Iterator begin() = 0;
  virtual Iterator end() = 0;
  virtual ConstIterator cbegin() const = 0;
  virtual ConstIterator cend() const = 0;

  virtual std::size_t size() const = 0;
};

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_QUADRATURE_SET_I_H_
