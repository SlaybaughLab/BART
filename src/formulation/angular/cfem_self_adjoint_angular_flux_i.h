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

  /*! \brief Initialize the formulation.
   * In general, this will pre-calculate matrix terms. The cell pointer is only
   * used to initialize the finite element object.
   * @param cell_ptr cell pointer for initialization.
   * @return Initialization token required for other calls
   */
  virtual InitializationToken Initialize(const formulation::CellPtr<dim>&) = 0;
};

} // namespace angular

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_ANGULAR_CFEM_SELF_ADJOINT_ANGULAR_FLUX_I_H_
