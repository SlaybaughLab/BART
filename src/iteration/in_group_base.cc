#include "in_group_base.h"
#include "../common/preconditioner_solver.h"

template <int dim>
InGroupBase<dim>::InGroupBase ()
:
IterationBase<dim> ()
{
  sol_ptr = build_linalg (prm);
}

template <int dim>
InGroupBase<dim>::~InGroupBase ()
{
}

template <int dim>
InGroupBase<dim>::solve_in_group ()
{
}

template class InGroupBase<2>;
template class InGroupBase<3>;
