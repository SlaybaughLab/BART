#ifndef BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_SET_H_
#define BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_SET_H_

#include "quadrature/angular/angular_quadrature_set_i.h"

namespace bart {

namespace quadrature {

namespace angular {

template <int dim>
class AngularQuadratureSet : public AngularQuadratureSetI<dim> {
 public:
  using typename AngularQuadratureSetI<dim>::AngleIndex;
  AngularQuadratureSet() = default;
  ~AngularQuadratureSet() = default;

  std::map<AngleIndex, QuadraturePoint<dim>> quadrature_points_map() const override {
    return quadrature_points_map_;
  }

  std::vector<QuadraturePoint<dim>> quadrature_points() const override {
    if (quadrature_points_.empty()) {
      quadrature_points_.resize(quadrature_points_map_.size());
      for (auto map_pair : quadrature_points_map_) {
        auto [angle_index, quadrature_point] = map_pair;
        quadrature_points_[angle_index] = quadrature_point;
      }
    }
    return quadrature_points_;
  }

  std::vector<Weight> quadrature_weights() const override {
    if (quadrature_weights_.empty()) {
      quadrature_weights_.resize(quadrature_points_map_.size());
      for (auto map_pair : quadrature_points_map_) {
        auto [angle_index, quadrature_point] = map_pair;
        quadrature_weights_[angle_index] = quadrature_point.first;
      }
    }
    return quadrature_weights_;
  }

  int total_quadrature_points() const override {
    return quadrature_points_map_.size();
  }

 protected:
  std::map<AngleIndex, QuadraturePoint<dim>> quadrature_points_map_ = {};

 private:
  mutable std::vector<QuadraturePoint<dim>> quadrature_points_ = {};
  mutable std::vector<Weight> quadrature_weights_ = {};
};

} // namespace angular

} // namespace quadrature

} // namespace bart

#endif // BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_SET_H_