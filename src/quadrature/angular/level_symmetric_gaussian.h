#ifndef BART_SRC_QUADRATURE_ANGULAR_LEVEL_SYMMETRIC_GAUSSIAN_H_
#define BART_SRC_QUADRATURE_ANGULAR_LEVEL_SYMMETRIC_GAUSSIAN_H_

#include "quadrature/angular/angular_quadrature_generator_i.h"

namespace bart {

namespace quadrature {

namespace angular {

class LevelSymmetricGaussian : AngularQuadratureGeneratorI<3> {
 public:
  explicit LevelSymmetricGaussian(quadrature::Order);
  std::vector<std::pair<CartesianPosition<3>, Weight>> GenerateSet() const;

  int order() const {
    return order_;
  }

 private:
  const int order_ = 0;
};


} // namespace angular

} // namespace quadrature

} //namespace bart

#endif //BART_SRC_QUADRATURE_ANGULAR_LEVEL_SYMMETRIC_GAUSSIAN_H_
