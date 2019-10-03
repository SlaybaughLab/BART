#include "quadrature/angular/level_symmetric_gaussian.h"

#include <deal.II/base/exceptions.h>

namespace bart {

namespace quadrature {

namespace angular {

LevelSymmetricGaussian::LevelSymmetricGaussian(bart::quadrature::Order order)
    : order_(order.get()) {
  AssertThrow(order_ >= 2,
              dealii::ExcMessage("Error in constructor of LevelSymmetricGaussian "
                                 "order must be >= 2"));
  AssertThrow(order_ % 2 == 0,
              dealii::ExcMessage("Error in constructor of LevelSymmetricGaussian "
                                 "order must be even"));
  AssertThrow(order_ <= 16,
              dealii::ExcMessage("Error in constructor of LevelSymmetricGaussian "
                                 "order must be <= 16"));
}

std::vector<std::pair<CartesianPosition<3>, Weight>> LevelSymmetricGaussian::GenerateSet() const {
  return std::vector<std::pair<CartesianPosition<3>, Weight>>();
}

} // namespace angular

} // namespace quadrature

} //namespace bart
