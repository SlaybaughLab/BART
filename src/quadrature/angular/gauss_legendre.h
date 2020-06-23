#ifndef BART_SRC_QUADRATURE_ANGULAR_GAUSS_LEGENDRE_H_
#define BART_SRC_QUADRATURE_ANGULAR_GAUSS_LEGENDRE_H_

#include "quadrature/quadrature_generator_i.h"

namespace bart {

namespace quadrature {

namespace angular {

/*! \brief Generates a 1D Gauss-Legendre quadrature set.
 *
 * The Gauss-Legendre polynomial set with $$n$$-points will integrate
 * polynomials up to degree $$2n-1$$ exactly. The points are generated on the
 * interval $$[0, 1]$$.
 *
 */
class GaussLegendre : public QuadratureGeneratorI<1> {
 public:
  //! Type for pairs of cartesian positions and weights
  using PositionWeightPairType = std::pair<CartesianPosition<1>, Weight>;

  /*! \brief Constructor.
   *
   * @param n_points number of points on $$[0, 1]$$
   */
  GaussLegendre(const int n_points);
  virtual ~GaussLegendre() = default;
  virtual std::vector<PositionWeightPairType> GenerateSet() const;
  //! \brief Returns the number of points
  virtual int order() const { return n_points_; };

 private:
  const int n_points_{0};
};

} // namespace angular

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_ANGULAR_GAUSS_LEGENDRE_H_
