#ifndef BART_SRC_FORMULATION_ASSEMBLE_SCALAR_H_
#define BART_SRC_FORMULATION_ASSEMBLE_SCALAR_H_

#include <memory>

#include "formulation/equation/diffusion.h"

#include "data/forward_declarations.h"
#include "data/system_scalar_fluxes.h"
#include "domain/definition.h"
#include "problem/parameter_types.h"
#include "utility/uncopyable.h"

namespace bart {

namespace formulation {

template <int dim>
class AssembleDiffusion : private utility::Uncopyable {
 public:
  using GroupNumber = int;
  using CellMatrix = dealii::FullMatrix<double>;
  using CellVector = dealii::Vector<double>;

  AssembleDiffusion(
      std::unique_ptr<equation::Diffusion<dim>> equation,
      std::unique_ptr<domain::Definition<dim>> domain,
      std::shared_ptr<data::ScalarSystemMatrixPtrs> system_matrix_ptrs,
      std::shared_ptr<data::ScalarRightHandSidePtrs> right_hand_side_ptrs,
      std::unique_ptr<data::ScalarRightHandSidePtrs> fixed_right_hand_side_ptrs,
      std::map<problem::Boundary, bool> reflective_boundary_map);
  ~AssembleDiffusion() = default;

//  void AssembleFixedBilinearTerms(GroupNumber group);
//  void AssembleFixedLinearTerms(GroupNumber group);
//  void AssembleVariableLinearTerms(GroupNumber group,
//                                   data::FluxVector &in_group_flux,
//                                   data::ScalarFluxPtrs &out_group_fluxes);

 private:
  // Unique pointers: equation and solver domain
  std::unique_ptr<equation::Diffusion<dim>> equation_;
  std::unique_ptr<domain::Definition<dim>> domain_;

  // Shared pointers -> System data to assemble into
  std::shared_ptr<data::ScalarSystemMatrixPtrs> system_matrix_ptrs_;
  std::shared_ptr<data::ScalarRightHandSidePtrs> right_hand_side_ptrs_;

  // Internal storage for fixed values
  std::unique_ptr<data::ScalarRightHandSidePtrs> fixed_right_hand_side_ptrs_;

  // Other
  std::map<problem::Boundary, bool> reflective_boundary_map_;
};

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_ASSEMBLE_SCALAR_H_