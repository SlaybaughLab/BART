#include "mg_base.h"

template <int dim>
MGBase<dim>::MGBase () : IterationBase<dim> ()
{
}

template <int dim>
MGBase<dim>::~MGBase ()
{
}

template <int dim>
void MGBase<dim>::do_iterations
(std::vector<PETScWrappers::MPI::SparseMatrix*> &sys_mats,
 std::vector<PETScWrappers::MPI::Vector*> &sys_rhses)
{
  this->initialize_equations ();
  multigroup_iterations ();
}

template <int dim>
void MGBase<dim>::multigroup_iterations
(std::vector<PETScWrappers::MPI::SparseMatrix*> &sys_mats,
 std::vector<PETScWrappers::MPI::Vector*> &sys_rhses)
{// this function needs to be overridden if JFNK is desired
  for (unsigned int g=0; g<n_group; ++g)
  {
    generate_group_rhses (sys_rhses, g);
    win_ptr->solve_in_group (sys_mats, g)
  }
}

template <int dim>
void MGBase<dim>::generate_group_rhses
(std::vector<PETScWrappers::MPI::Vector*> &group_rhses, unsigned int &g)
{// invoke win_ptr to assemble rhs per specific MG solver in derived class
}

/*
template <int dim>
void MGBase<dim>::iterate_over_groups
(std::vector<PETScWrappers::MPI::Vector*> &group_rhses)
{
  for (unsigned int i=0; i<n_group; ++g)
  {
    generate_group_rhses (group_rhses, g);
    win_ptr->solve_in_group (g);
  }
}
*/

template class MGBase<2>;
template class MGBase<3>;
