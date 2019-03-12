#include "formulation/assemble_scalar.h"

namespace bart {

namespace formulation {

template <int dim>
AssembleScalar<dim>::AssembleScalar(
    std::unique_ptr<equation::TransportScalar<dim>> equation,
    std::unique_ptr<domain::Definition<dim>> domain,
    std::shared_ptr<data::SystemScalarFluxes> scalar_fluxes,
    std::shared_ptr<data::ScalarSystemMatrices> system_matrices,
    std::shared_ptr<data::RightHandSideVector> right_hand_side)
    : equation_(std::move(equation)),
      domain_(std::move(domain)),
      scalar_fluxes_(scalar_fluxes),
      system_matrices_(system_matrices),
      right_hand_side_(right_hand_side) {


    }

template class AssembleScalar<1>;
template class AssembleScalar<2>;
template class AssembleScalar<3>;

} // namespace formulation

} // namespace bart

