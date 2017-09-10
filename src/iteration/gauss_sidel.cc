#include "mg_base.h"

template <int dim>
GaussSidel<dim>::GaussSidel (ParameterHandler &prm)
:
MGBase<dim> (prm)
{
}

template <int dim>
GaussSidel<dim>::~GaussSidel ()
{
}

template <int dim>
void GaussSidel<dim>::mg_iterations
(std::vector<Vector<double> > &sflxes_proc,
 std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs)
{
  // TODO: fill this up
}

template class GaussSidel<2>;
template class GaussSidel<3>;
