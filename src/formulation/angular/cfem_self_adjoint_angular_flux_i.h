#ifndef BART_SRC_FORMULATION_ANGULAR_CFEM_SELF_ADJOINT_ANGULAR_FLUX_I_H_
#define BART_SRC_FORMULATION_ANGULAR_CFEM_SELF_ADJOINT_ANGULAR_FLUX_I_H_

#include <deal.II/lac/full_matrix.h>
#include <deal.II/dofs/dof_accessor.h>

#include "formulation/formulation_types.h"

namespace bart {

namespace formulation {

namespace angular {

template <int dim>
class CFEMSelfAdjointAngularFluxI {
 public:
  struct InitializationToken{};

  virtual ~CFEMSelfAdjointAngularFluxI() = default;
};

} // namespace angular

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_ANGULAR_CFEM_SELF_ADJOINT_ANGULAR_FLUX_I_H_
