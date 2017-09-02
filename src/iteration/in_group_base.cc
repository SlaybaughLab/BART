#include "in_group_base.h"
#include "../common/preconditioner_solver.h"

template <int dim>
InGroupBase<dim>::InGroupBase ()
:
IterationBase<dim> (),
err_phi_tol(1.0e-6)
{
}

template <int dim>
InGroupBase<dim>::~InGroupBase ()
{
}

template <int dim>
InGroupBase<dim>::solve_in_group
(std_cxx11::shared_ptr<EquationBase<dim> > equ_ptr,
 unsigned int &g)
{
}

template <int dim>
SourceIteration<dim>::SourceIteration
:
InGroupBase<dim> ()
{
}

template <int dim>
SourceIteration<dim>::solve_in_group
(std_cxx11::shared_ptr<EquationBase<dim> > equ_ptr,
 unsigned int &g)
{
  double err = 1.0;
  while (err<this->err_phi_tol)
}

template class InGroupBase<2>;
template class InGroupBase<3>;
