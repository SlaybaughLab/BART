#include "in_group_base.h"

template <int dim>
InGroupBase<dim>::InGroupBase ()
:
IterationBase<dim> (),
err_phi_tol(1.0e-6)
{
  sflx_proc_old.resize (1);
}

template <int dim>
InGroupBase<dim>::~InGroupBase ()
{
}

template <int dim>
InGroupBase<dim>::solve_in_group
(std::vector<Vector<double> > &sflx_proc,
 std_cxx11::shared_ptr<EquationBase<dim> > equ_ptrs,
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
(std::vector<Vector<double> > &sflx_proc,
 std_cxx11::shared_ptr<EquationBase<dim> > equ_ptrs,
 unsigned int &g)
{
  double err = 1.0;
  while (err>this->err_phi_tol)
  {
    // generate rhs for group g
    equ_ptrs[0]->assemble_linear_form (sflx_proc, g);
    // solve all the directions in group g
    equ_ptrs[0]->solve_in_group (g);
    // generate moments
    equ_ptrs[0]->generate (sflx_proc[g], sflx_proc_old, g);
    // calculate the difference of moments for convergence check
    err = this->estimate_phi_diff (sflx_proc[g], sflx_proc_old);
  }
}

template class InGroupBase<2>;
template class InGroupBase<3>;
