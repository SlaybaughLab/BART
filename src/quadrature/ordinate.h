#ifndef BART_SRC_QUADRATURE_ORDINATE_H_
#define BART_SRC_QUADRATURE_ORDINATE_H_

#include "quadrature/ordinate_i.h"

namespace bart {

namespace quadrature {

template <int dim>
class Ordinate : public OrdinateI<dim> {
 public:
  virtual ~Ordinate() = default;
};

} // namespace quadrature

} //namespace bart

#endif //BART_SRC_QUADRATURE_ORDINATE_H_
