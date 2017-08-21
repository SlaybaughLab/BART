#include "source_iteration.h"

SourceIteration<dim>::SourceIteration ()
:
IterationBase<dim> ()
{
}

SourceIteration<dim>::~SourceIteration ()
{
}

SourceIteration<dim>::solve_in_group ()
{
  // complete source iteration
  for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
  {
    alg_ptr->ho_solve ();
  }
}

template class SourceIteration<2>;
template class SourceIteration<3>;
