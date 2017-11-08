#include <deal.II/base/index_set.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/cell_id.h>

#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/solver_bicgstab.h>

#include <algorithm>

#include "iterations.h"

using namespace dealii;

template <int dim>
Iterations<dim>::Iterations
(const ParameterHandler &prm)
:
transport_name(prm.get("transport model")),
is_eigen_problem(prm.get_bool("do eigenvalue calculations"))
{
}

template <int dim>
Iterations<dim>::~Iterations ()
{
}

template <int dim>
void Iterations<dim>::solve_problems
(std::vector<Vector<double> > &sflxes_proc,
 std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs,
 std_cxx11::shared_ptr<IGBase<dim> > ig_ptr,
 std_cxx11::shared_ptr<MGBase<dim> > mg_ptr,
 std_cxx11::shared_ptr<EigenBase<dim> > eig_ptr)
{
  if (is_eigen_problem)
  {
    eig_ptr->do_iterations (sflxes_proc, equ_ptrs, ig_ptr, mg_ptr);
    keff = eig_ptr->get_keff ();
  }
  else
    mg_ptr->do_iterations (sflxes_proc, equ_ptrs, ig_ptr);
}

// explicit instantiation to avoid linking error
template class Iterations<2>;
template class Iterations<3>;
