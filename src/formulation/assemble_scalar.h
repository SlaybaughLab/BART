#ifndef BART_SRC_FORMULATION_ASSEMBLE_SCALAR_H_
#define BART_SRC_FORMULATION_ASSEMBLE_SCALAR_H_

#include <memory>

#include "data/forward_declarations.h"
#include "data/system_scalar_fluxes.h"
#include "domain/definition.h"
#include "formulation/equation/transport_scalar.h"
#include "problem/parameter_types.h"
#include "utility/uncopyable.h"

namespace bart {

namespace formulation {

template <int dim>
class AssembleScalar : private utility::Uncopyable {
 public:
  using GroupNumber = int;
  using CellMatrix = dealii::FullMatrix<double>;
  using CellVector = dealii::Vector<double>;

  AssembleScalar(
      std::unique_ptr<equation::TransportScalar<dim>> equation,
      std::unique_ptr<domain::Definition<dim>> domain,
      std::shared_ptr<data::SystemScalarFluxes> scalar_fluxes,
      std::shared_ptr<data::ScalarSystemMatrices> system_matrices,
      std::shared_ptr<data::ScalarRightHandSideVectors> right_hand_side,
      std::unique_ptr<data::ScalarRightHandSideVectors> fixed_right_hand_side,
      std::map<problem::Boundary, bool> reflective_boundary_map);
  ~AssembleScalar() = default;

  void AssembleBilinearTerms(GroupNumber group);
  void AssembleFixedSourceLinearTerm(GroupNumber group);

 private:
  // Unique pointers: equation and solver domain
  std::unique_ptr<equation::TransportScalar<dim>> equation_;
  std::unique_ptr<domain::Definition<dim>> domain_;

  // Shared pointers -> System data to assemble into
  std::shared_ptr<data::SystemScalarFluxes> scalar_fluxes_;
  std::shared_ptr<data::ScalarSystemMatrices> system_matrices_;
  std::shared_ptr<data::ScalarRightHandSideVectors> right_hand_side_;

  // Internal storage for fixed RHS values
  std::unique_ptr<data::ScalarRightHandSideVectors> fixed_right_hand_side_;

  // Other
  std::map<problem::Boundary, bool> reflective_boundary_map_;
};

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_ASSEMBLE_SCALAR_H_