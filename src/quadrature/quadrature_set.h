#ifndef BART_SRC_QUADRATURE_QUADRATURE_SET_H_
#define BART_SRC_QUADRATURE_QUADRATURE_SET_H_

#include "quadrature/quadrature_set_i.h"

namespace bart {

namespace quadrature {

/*! \brief Default implementation of the QuadratureSet.
 *
 * @tparam dim spatial dimension
 */
template <int dim>
class QuadratureSet : public QuadratureSetI<dim> {
 public:
  virtual ~QuadratureSet() = default;

  bool AddPoint(std::shared_ptr<QuadraturePointI<dim>>);
  size_t size() const override {
    return quadrature_point_ptrs_.size();
  };

 protected:
  std::set<std::shared_ptr<QuadraturePointI<dim>>> quadrature_point_ptrs_;
};

} // namespace quadrature

} //namespace bart

#endif //BART_SRC_QUADRATURE_QUADRATURE_SET_H_
