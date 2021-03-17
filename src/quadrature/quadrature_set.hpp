#ifndef BART_SRC_QUADRATURE_QUADRATURE_SET_HPP_
#define BART_SRC_QUADRATURE_QUADRATURE_SET_HPP_

#include "quadrature/quadrature_set_i.hpp"
#include "quadrature/utility/quadrature_utilities.h"

namespace bart::quadrature {

/*! \brief Default implementation of the QuadratureSet.
 *
 * @tparam dim spatial dimension
 */
template <int dim>
class QuadratureSet : public QuadratureSetI<dim> {
 public:
  using QuadraturePoint = QuadraturePointI<dim>;
  using typename QuadratureSetI<dim>::Iterator, typename QuadratureSetI<dim>::ConstIterator;
  virtual ~QuadratureSet() = default;

  auto AddPoint(std::shared_ptr<QuadraturePoint>) -> bool;
  [[nodiscard]] auto GetBoundaryReflection(const std::shared_ptr<QuadraturePoint>& quadrature_point_to_reflect,
                                           problem::Boundary) const -> std::shared_ptr<QuadraturePoint> override;
  auto SetReflection(std::shared_ptr<QuadraturePoint>, std::shared_ptr<QuadraturePoint>) -> void;

  [[nodiscard]] auto GetReflection(std::shared_ptr<QuadraturePoint>) const -> std::shared_ptr<QuadraturePoint> override;
  [[nodiscard]] auto GetReflectionIndex(std::shared_ptr<QuadraturePoint>) const -> std::optional<int> override;
  [[nodiscard]] auto GetQuadraturePoint(QuadraturePointIndex) const -> std::shared_ptr<QuadraturePoint> override;
  [[nodiscard]] auto GetQuadraturePointIndex(std::shared_ptr<QuadraturePoint>) const -> int override;
  [[nodiscard]] auto quadrature_point_indices() const -> std::set<int> override { return quadrature_point_indices_; };

  auto begin() -> Iterator override { return quadrature_point_ptrs_.begin(); };
  auto end() -> Iterator override { return quadrature_point_ptrs_.end(); };
  auto cbegin() const -> ConstIterator override { return quadrature_point_ptrs_.cbegin(); };
  auto cend() const -> ConstIterator override { return quadrature_point_ptrs_.cend(); };
  [[nodiscard]] auto size() const -> std::size_t override {return quadrature_point_ptrs_.size(); };

 protected:
  //! Set of quadrature point indices
  std::set<int> quadrature_point_indices_{};
  //! Set of actual quadrature points
  std::set<std::shared_ptr<QuadraturePoint>, utility::quadrature_point_compare<dim>> quadrature_point_ptrs_;
  //! Mapping of points to their reflections
  std::map<std::shared_ptr<QuadraturePoint>, std::shared_ptr<QuadraturePoint>> reflection_map_;
  //! Mapping of quadrature point indices to quadrature points
  std::map<int, std::shared_ptr<QuadraturePoint>> index_to_quadrature_point_map_{};
  //! Mapping of indices to quadrature point
  std::map<std::shared_ptr<QuadraturePoint>, int> quadrature_point_to_index_map_{};
  //! Internal mapping of boundary reflections
  mutable std::map<problem::Boundary, std::map<std::shared_ptr<QuadraturePoint>, std::shared_ptr<QuadraturePoint>>>
      boundary_reflections_;
};

} // namespace bart::quadrature

#endif //BART_SRC_QUADRATURE_QUADRATURE_SET_HPP_
