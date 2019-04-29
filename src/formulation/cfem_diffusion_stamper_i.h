#ifndef BART_SRC_FORMULATION_CFEM_DIFFUSION_STAMPER_I_H_
#define BART_SRC_FORMULATION_CFEM_DIFFUSION_STAMPER_I_H_

#include "formulation/stamper_i.h"

namespace bart {

namespace formulation {

template <int dim>
class CFEM_DiffusionStamperI : public StamperI<dim> {
 public:
  virtual ~CFEM_DiffusionStamperI() = default;
};

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_CFEM_DIFFUSION_STAMPER_I_H_