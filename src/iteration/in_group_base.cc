#include "in_group_base.h"

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
(std::vector<Vector<double> > &sflxes_proc,
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
(std::vector<Vector<double> > &sflxes_proc,
 std_cxx11::shared_ptr<EquationBase<dim> > equ_ptrs,
 unsigned int &g)
{
  double err = 1.0;
  while (err>this->err_phi_tol)
  {
    // generate rhs for group g
    equ_ptrs[0]->assemble_linear_form (sflxes_proc, g);
    // solve all the directions in group g
    equ_ptrs[0]->solve_in_group (g);
    // generate moments
    equ_ptrs[0]->generate (sflxes_proc[g], this->sflx_proc_prev_ig, g);
    // calculate the difference of moments for convergence check
    err = this->estimate_phi_diff (sflxes_proc[g], this->sflx_proc_prev_ig);
  }
}

template class InGroupBase<2>;
template class InGroupBase<3>;
