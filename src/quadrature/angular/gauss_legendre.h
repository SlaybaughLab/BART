#ifndef BART_SRC_QUADRATURE_ANGULAR_GAUSS_LEGENDRE_H_
#define BART_SRC_QUADRATURE_ANGULAR_GAUSS_LEGENDRE_H_

#include "quadrature/quadrature_generator_i.h"

namespace bart {

namespace quadrature {

namespace angular {

/*! \brief 1D Gauss-Legendre quadrature set.
 *
 */
class GaussLegendre : public QuadratureGeneratorI<1> {
 public:
  using PositionWeightPairType = std::pair<CartesianPosition<1>, Weight>;

  GaussLegendre(const int n_points) : n_points_(n_points) {}
  virtual ~GaussLegendre() = default;
  virtual std::vector<PositionWeightPairType> GenerateSet() const;
  virtual int order() const { return n_points_; };

 private:
  const int n_points_ = 0;
};

} // namespace angular

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_ANGULAR_GAUSS_LEGENDRE_H_
