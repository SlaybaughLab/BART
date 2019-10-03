#include "quadrature/angular/level_symmetric_gaussian.h"

namespace bart {

namespace quadrature {

namespace angular {

LevelSymmetricGaussian::LevelSymmetricGaussian(bart::quadrature::Order order)
    : order_(order.get()) {}

std::vector<std::pair<CartesianPosition<3>, Weight>> LevelSymmetricGaussian::GenerateSet() const {
  return std::vector<std::pair<CartesianPosition<3>, Weight>>();
}

} // namespace angular

} // namespace quadrature

} //namespace bart
