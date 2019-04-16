#ifndef BART_SRC_FORMULATION_SCALAR_CFEM_DIFFUSION_I_H_
#define BART_SRC_FORMULATION_SCALAR_CFEM_DIFFUSION_I_H_

#include "formulation/scalar/cfem_i.h"

namespace bart {

namespace formulation {

namespace scalar {

template <int dim>
class CFEM_DiffusionI : public CFEM_I {
 public:
  virtual ~CFEM_DiffusionI() = default;
};

} // namespace scalar

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_SCALAR_CFEM_DIFFUSION_I_H_