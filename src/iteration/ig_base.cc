#include "ig_base.h"

template <int dim>
IGBase<dim>::IGBase (const ParameterHandler &prm)
:
IterationBase<dim> (prm),
err_phi_tol(1.0e-6)
{
}

template <int dim>
IGBase<dim>::~IGBase ()
{
}

template <int dim>
void IGBase<dim>::solve_in_group
(std::vector<Vector<double> > &sflxes_proc,
 std_cxx11::shared_ptr<EquationBase<dim> > equ_ptr,
 unsigned int &g)
{
  // the default is for diffusion like system, SPN and PN
  equ_ptr->solve_in_group (g);
}

template <int dim>
SourceIteration<dim>::SourceIteration (const ParameterHandler &prm)
:
IGBase<dim> (prm)
{
}

template <int dim>
SourceIteration<dim>::~SourceIteration ()
{
}

template <int dim>
void SourceIteration<dim>::solve_in_group
(std::vector<Vector<double> > &sflxes_proc,
 std_cxx11::shared_ptr<EquationBase<dim> > equ_ptr,
 unsigned int &g)
{
  double err = 1.0;
  int iter = 0;
  while (err>this->err_phi_tol)
  {
    // generate rhs for group g
    equ_ptr->assemble_linear_form (sflxes_proc, g);
    // solve all the directions in group g
    equ_ptr->solve_in_group (g);
    // generate moments
    equ_ptr->generate_moments (sflxes_proc[g], this->sflx_proc_prev_ig, g);
    // calculate the difference of moments for convergence check
    err = this->estimate_phi_diff (sflxes_proc[g], this->sflx_proc_prev_ig);
    this->pcout << "SI error: " << err << ", Iter: " << iter++ << std::endl;
  }
}

template class IGBase<2>;
template class IGBase<3>;
template class SourceIteration<2>;
template class SourceIteration<3>;
