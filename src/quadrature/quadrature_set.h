#ifndef BART_SRC_QUADRATURE_QUADRATURE_SET_H_
#define BART_SRC_QUADRATURE_QUADRATURE_SET_H_

#include "quadrature/quadrature_set_i.h"
#include "quadrature/utility/quadrature_utilities.h"

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

  void SetReflection(std::shared_ptr<QuadraturePointI<dim>>,
                     std::shared_ptr<QuadraturePointI<dim>>);

  std::shared_ptr<QuadraturePointI<dim>> GetReflection(
      std::shared_ptr<QuadraturePointI<dim>>) const override;

  typename QuadratureSetI<dim>::Iterator begin() override {
    return quadrature_point_ptrs_.begin();
  };
  typename QuadratureSetI<dim>::Iterator end() override {
    return quadrature_point_ptrs_.end();
  }

  typename QuadratureSetI<dim>::ConstIterator cbegin() const override {
    return quadrature_point_ptrs_.cbegin();
  };

  typename QuadratureSetI<dim>::ConstIterator cend() const override {
    return quadrature_point_ptrs_.cend();
  };

  size_t size() const override {
    return quadrature_point_ptrs_.size();
  };

 protected:
  std::set<std::shared_ptr<QuadraturePointI<dim>>,
           utility::quadrature_point_compare<dim>> quadrature_point_ptrs_;
  //std::set<std::shared_ptr<QuadraturePointI<dim>>> quadrature_point_ptrs_;
  std::map<std::shared_ptr<QuadraturePointI<dim>>,
           std::shared_ptr<QuadraturePointI<dim>>> reflection_map_;
};

} // namespace quadrature

} //namespace bart

#endif //BART_SRC_QUADRATURE_QUADRATURE_SET_H_
