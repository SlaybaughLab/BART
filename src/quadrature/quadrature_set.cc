#include "quadrature/quadrature_set.h"

#include <algorithm>

namespace bart {

namespace quadrature {

template<int dim>
bool QuadratureSet<dim>::AddPoint(
    std::shared_ptr<QuadraturePointI<dim>> new_point_ptr) {
  AssertThrow(new_point_ptr != nullptr,
      dealii::ExcMessage("Error in AddPoint, pointer is null"))
  auto status = quadrature_point_ptrs_.insert(new_point_ptr);

  auto max_index = std::max_element(quadrature_point_indices_.begin(),
                                    quadrature_point_indices_.end());

  if (status.second) {
    int new_index = 0;

    if (max_index != quadrature_point_indices_.end())
      new_index = *max_index + 1;

    // Update sets and maps
    quadrature_point_indices_.insert(new_index);
    index_to_quadrature_point_map_.insert_or_assign(new_index, new_point_ptr);
    quadrature_point_to_index_map_.insert_or_assign(new_point_ptr, new_index);
  }

  return status.second;
}

template<int dim>
std::shared_ptr<QuadraturePointI<dim>> QuadratureSet<dim>::GetBoundaryReflection(
    const std::shared_ptr<QuadraturePointI<dim>> &quadrature_point_to_reflect,
    problem::Boundary boundary) const {
  using Boundary = problem::Boundary;

  std::shared_ptr<QuadraturePointI<dim>> return_ptr = nullptr;

  if (auto boundary_mapping = boundary_reflections_.find(boundary);
      boundary_mapping != boundary_reflections_.end()) {
    if (auto reflection_pair = boundary_mapping->second.find(quadrature_point_to_reflect);
        reflection_pair != boundary_mapping->second.end()) {
      return_ptr = reflection_pair->second;
    }
  } else {
    boundary_reflections_.insert({boundary, {}});
  }

  if (return_ptr == nullptr) {
    auto position = quadrature_point_to_reflect->cartesian_position();

    switch (boundary) {
      case Boundary::kXMax: case Boundary::kXMin: {
        position.at(0) *= -1;
        break;
      }
      case Boundary::kYMin: case Boundary::kYMax: {
        position.at(1) *= -1;
        break;
      }
      case Boundary::kZMin: case Boundary::kZMax: {
        position.at(2) *= -1;
        break;
      }
    }

    for (const auto& quadrature_point : this->quadrature_point_ptrs_) {
      if (quadrature_point->cartesian_position() == position) {
        return_ptr = quadrature_point;
        break;
      }
    }

    boundary_reflections_.at(boundary).insert({quadrature_point_to_reflect,
                                               return_ptr});
    boundary_reflections_.at(boundary).insert({return_ptr,
                                               quadrature_point_to_reflect});

    AssertThrow(return_ptr != nullptr,
                dealii::ExcMessage("GetBoundaryReflection returned null reflection"))
  }

  return return_ptr;
}

template<int dim>
void QuadratureSet<dim>::SetReflection(
    std::shared_ptr<QuadraturePointI<dim>> first_point,
    std::shared_ptr<QuadraturePointI<dim>> second_point) {

  AssertThrow(first_point != second_point,
              dealii::ExcMessage("Error in SetReflection: both points are the "
                                 "same"));

  bool both_points_in_set = (quadrature_point_ptrs_.count(first_point) == 1) &&
      (quadrature_point_ptrs_.count(second_point) == 1);

  AssertThrow(both_points_in_set,
      dealii::ExcMessage("Error in SetReflection: one or both points are not "
                         "in the quadrature set"));

  for (auto& point : {first_point, second_point}) {
    try {
      // If the first point is already present, we will delete the entry and
      // the entry for its previous reflection
      auto old_reflection = reflection_map_.at(point);
      reflection_map_.erase(old_reflection);
    } catch (const std::out_of_range&) {}
  }

  reflection_map_.insert_or_assign(first_point, second_point);
  reflection_map_.insert_or_assign(second_point, first_point);
}

template<int dim>
std::shared_ptr<QuadraturePointI<dim>> QuadratureSet<dim>::GetReflection(
    std::shared_ptr<QuadraturePointI<dim>> quadrature_point) const {
  try {
    return reflection_map_.at(quadrature_point);
  } catch (const std::out_of_range&) {
    AssertThrow(quadrature_point_ptrs_.count(quadrature_point) == 1,
                dealii::ExcMessage("Error in GetReflection, quadrature point "
                                   "is not in quadrature set."))
    return nullptr;
  }
}

template<int dim>
std::optional<int> QuadratureSet<dim>::GetReflectionIndex(
    std::shared_ptr<QuadraturePointI<dim>> quadrature_point_ptr) const {

  std::optional<int> reflection_index;

  auto reflection_ptr = GetReflection(quadrature_point_ptr);

  if (reflection_ptr != nullptr)
    reflection_index = quadrature_point_to_index_map_.at(reflection_ptr);

  return reflection_index;
}

template<int dim>
std::shared_ptr<QuadraturePointI<dim>> QuadratureSet<dim>::GetQuadraturePoint(
    QuadraturePointIndex index) const {
  return index_to_quadrature_point_map_.at(index.get());
}
template<int dim>
int QuadratureSet<dim>::GetQuadraturePointIndex(
    std::shared_ptr<QuadraturePointI<dim>> quadrature_point) const {
  return quadrature_point_to_index_map_.at(quadrature_point);
}


template class QuadratureSet<1>;
template class QuadratureSet<2>;
template class QuadratureSet<3>;

} // namespace quadrature

} // namespace bart

