#include "self_adjoint_angular_flux.h"

template <int dim>
SelfAdjointAngularFlux<dim>::SelfAdjointAngularFlux(
    const std::string equation_name,
    const dealii::ParameterHandler &prm,
    std::shared_ptr<FundamentalData<dim>> &data_ptr)
    : EquationBase<dim>(equation_name, prm, data_ptr) {}

// cfem

template <int dim>
void SelfAdjointAngularFlux<dim>::PreassembleCellMatrices () {}
