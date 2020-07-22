#include "quadrature/angular/gauss_legendre.h"

#include <deal.II/base/quadrature_lib.h>

namespace bart {

namespace quadrature {

namespace angular {

GaussLegendre::GaussLegendre(const int n_points)
    : n_points_(n_points) {
  AssertThrow(n_points > 0,
      dealii::ExcMessage("Error in constructor of GaussLegendre, n_points must "
                         "be greater than or equal to 0"))
  this->set_description("Gauss-Legendre quadrature (1D)");
}


std::vector<GaussLegendre::PositionWeightPairType>
    GaussLegendre::GenerateSet() const {
  dealii::QGauss<1> gaussian_quadrature(n_points_);
  std::vector<PositionWeightPairType> return_vector;

  auto points = gaussian_quadrature.get_points();
  auto weights = gaussian_quadrature.get_weights();

  for (int i = 0; i < static_cast<int>(gaussian_quadrature.size()); ++i) {
    CartesianPosition<1> x_position({points.at(i)[0]});
    Weight weight(2*M_PI*weights.at(i));
    return_vector.push_back({x_position, weight});
  }

  return return_vector;
}

} // namespace angular

} // namespace quadrature

} // namespace bart