#include "equation_base.h"

template <int dim>
EquationBase<dim>::EquationBase (dealii::ParameterHandler &prm,
    std::unique_ptr<FundamentalData<dim>> &fund_dat_ptr)
    : dat_ptr(fund_dat_ptr) {
}

template <int dim>
EquationBase<dim>::~EquationBase () {}
