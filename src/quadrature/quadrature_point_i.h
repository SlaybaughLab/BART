#ifndef BART_SRC_QUADRATURE_QUADRATURE_POINT_I_H_
#define BART_SRC_QUADRATURE_QUADRATURE_POINT_I_H_

namespace bart {

namespace quadrature {

template <int dim>
class QuadraturePointI {
 public:
  virtual ~QuadraturePointI() = default;
};

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_QUADRATURE_POINT_I_H_
