#include "quadrature/utility/quadrature_utilities.h"
#include "quadrature/quadrature_set_i.h"

#include <functional>

namespace bart {

namespace quadrature {

namespace utility {

template <int dim>
std::array<double, dim> ReflectAcrossOrigin(const OrdinateI<dim>& ordinate) {

  auto position = ordinate.cartesian_position();

  std::transform(position.begin(), position.end(),
                 position.begin(), [](double val){return -val;});

  return position;
}

template <>
std::vector<std::pair<CartesianPosition<1>, Weight>> GenerateAllPositiveX<1>(
    const std::vector<std::pair<CartesianPosition<1>, Weight>>& to_distribute) {
  return to_distribute;
}

template <>
std::vector<std::pair<CartesianPosition<2>, Weight>> GenerateAllPositiveX<2>(
    const std::vector<std::pair<CartesianPosition<2>, Weight>>& to_distribute) {
  auto quadrature_pairs = to_distribute;
  for (auto [position, weight] : to_distribute) {
    auto& y = position.get().at(1);
    if (y != 0) {
      y *= -1;
      quadrature_pairs.emplace_back(CartesianPosition<2>(position),
                                    Weight(weight));
    }
  }
  return quadrature_pairs;
}

template <>
std::vector<std::pair<CartesianPosition<3>, Weight>> GenerateAllPositiveX<3>(
    const std::vector<std::pair<CartesianPosition<3>, Weight>>& to_distribute) {
  auto quadrature_pairs = to_distribute;
  for (auto [position, weight] : to_distribute) {
    auto& y = position.get().at(1);
    auto& z = position.get().at(2);
    if (y != 0) {
      y *= -1;
      quadrature_pairs.emplace_back(CartesianPosition<3>(position),
                                    Weight(weight));
    }
    if (z != 0) {
      z *= -1;
      quadrature_pairs.emplace_back(CartesianPosition<3>(position),
                                    Weight(weight));
      if (y != 0) {
        y *= -1;
        quadrature_pairs.emplace_back(CartesianPosition<3>(position),
                                      Weight(weight));
      }
    }
  }
  return quadrature_pairs;
}

template std::array<double, 1> ReflectAcrossOrigin<1>(const OrdinateI<1>&);
template std::array<double, 2> ReflectAcrossOrigin<2>(const OrdinateI<2>&);
template std::array<double, 3> ReflectAcrossOrigin<3>(const OrdinateI<3>&);

} // namespace utility

} // namespace quadrature

} // namespace bart