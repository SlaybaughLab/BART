#include "quadrature/angular/level_symmetric_gaussian.h"

#include <deal.II/base/exceptions.h>
#include <deal.II/base/quadrature_lib.h>

namespace bart {

namespace quadrature {

namespace angular {

LevelSymmetricGaussian::LevelSymmetricGaussian(bart::quadrature::Order order)
    : order_(order.get()) {
  AssertThrow(order_ >= 2,
              dealii::ExcMessage("Error in constructor of "
                                 "LevelSymmetricGaussian order must be >= 2"))
  AssertThrow(order_ % 2 == 0,
              dealii::ExcMessage("Error in constructor of "
                                 "LevelSymmetricGaussian order must be even"))
  AssertThrow(order_ <= 16,
              dealii::ExcMessage("Error in constructor of "
                                 "LevelSymmetricGaussian order must be <= 16"))
  this->set_description("Level-symmetric-type Gaussian quadrature (3D)");
}

std::vector<std::pair<CartesianPosition<3>, Weight>>
LevelSymmetricGaussian::GenerateSet() const {
  std::vector<std::pair<CartesianPosition<3>, Weight>> generated_set;
  const int n_points = order_/2;

  // Gaussian quadrature for theta, the lowest level has order_ points but we
  // will only use half of the quadrature
  dealii::QGauss<1> gaussian_quadrature(order_);

  // For each quadrature point in the lowest level, we will create equally
  // spaced points in phi
  for (int level = 0; level < n_points; ++level) {
    int n_points_this_level = level + 1;
    // mu is just the gaussian quadrature point
    double mu = 1 - gaussian_quadrature.point(level)[0] * 2;
    // points in phi will be equally spaced from 0 to pi/2
    double dphi = M_PI/(2*n_points_this_level);
    // weights on this level are equal to the gaussian times 2PI/total points in
    // 2PI. Here we have n_points_this_level*4 total points
    double weight =
        gaussian_quadrature.weight(level) * M_PI/(n_points_this_level);

    for (int j = 0; j < n_points_this_level; ++j) {
      double phi = (j + 0.5)*dphi;
      double x = std::sqrt(1 - mu*mu) * std::cos(phi);
      double y = std::sqrt(1 - mu*mu) * std::sin(phi);
      double z = mu;
      std::array<double, 3> position{x,y,z};
      generated_set.emplace_back(std::make_pair(CartesianPosition<3>(position),
                                                Weight(weight)));
    }
  }
  return generated_set;
}

} // namespace angular

} // namespace quadrature

} //namespace bart