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
  std::shared_ptr<QuadraturePointI<dim>> GetBoundaryReflection(
      const std::shared_ptr<QuadraturePointI<dim>>& quadrature_point_to_reflect,
      problem::Boundary boundary) const override;
  std::shared_ptr<QuadraturePointI<dim>> GetReflection(
      std::shared_ptr<QuadraturePointI<dim>>) const override;
  std::optional<int> GetReflectionIndex(
      std::shared_ptr<QuadraturePointI<dim>>) const override;
  std::shared_ptr<QuadraturePointI<dim>> GetQuadraturePoint(
      QuadraturePointIndex index) const override;
  int GetQuadraturePointIndex(
      std::shared_ptr<QuadraturePointI<dim>> quadrature_point) const override;
  std::set<int> quadrature_point_indices() const override {
    return quadrature_point_indices_; };
  typename QuadratureSetI<dim>::Iterator begin() override {
    return quadrature_point_ptrs_.begin(); };
  typename QuadratureSetI<dim>::Iterator end() override {
    return quadrature_point_ptrs_.end(); };
  typename QuadratureSetI<dim>::ConstIterator cbegin() const override {
    return quadrature_point_ptrs_.cbegin(); };
  typename QuadratureSetI<dim>::ConstIterator cend() const override {
    return quadrature_point_ptrs_.cend(); };
  size_t size() const override {return quadrature_point_ptrs_.size(); };

 protected:
  //! Set of quadrature point indices
  std::set<int> quadrature_point_indices_ = {};
  //! Set of actual quadrature points
  std::set<std::shared_ptr<QuadraturePointI<dim>>,
           utility::quadrature_point_compare<dim>> quadrature_point_ptrs_;
  //! Mapping of points to their reflections
  std::map<std::shared_ptr<QuadraturePointI<dim>>,
           std::shared_ptr<QuadraturePointI<dim>>> reflection_map_;
  //! Mapping of quadrature point indices to quadrature points
  std::map<int, std::shared_ptr<QuadraturePointI<dim>>>
      index_to_quadrature_point_map_ = {};
  //! Mapping of indices to quadrature point
  std::map<std::shared_ptr<QuadraturePointI<dim>>, int>
      quadrature_point_to_index_map_ = {};
};

} // namespace quadrature

} //namespace bart

#endif //BART_SRC_QUADRATURE_QUADRATURE_SET_H_
