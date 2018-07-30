#include "ig_base.h"

template <int dim>
IGBase<dim>::IGBase (const dealii::ParameterHandler &prm,
    std::shared_ptr<FundamentalData<dim>> &dat_ptr)
    :
    IterationBase<dim> (prm, dat_ptr),
    err_phi_tol_(1.0e-6) {}

template <int dim>
IGBase<dim>::~IGBase () {}

template <int dim>
void IGBase<dim>::IGIterations (std::unique_ptr<EquationBase<dim>> &equ_ptr,
    const int &g) {
  // the default is for diffusion like system, SPN and PN
  equ_ptr->SolveInGroup (g);
}

template class IGBase<1>;
template class IGBase<2>;
template class IGBase<3>;
