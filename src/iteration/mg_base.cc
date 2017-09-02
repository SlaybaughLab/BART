#include "mg_base.h"

template <int dim>
MGBase<dim>::MGBase ()
: IterationBase<dim> (),
err_phi_tol(1.0e-5)
{
}

template <int dim>
MGBase<dim>::~MGBase ()
{
}

template <int dim>
void MGBase<dim>::do_iterations
(std::vector<typename DoFHandler<dim>::active_cell_iterator> &local_cells,
 std::vector<bool> &is_cell_at_bd,
 std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs)
{
  // assemble bilinear form of transport equation
  equ_ptrs[0]->assemble_bilinear_forms (local_cells, is_cell_at_bd);
  // assemble NDA bilinear form if do_nda
  if (do_nda)
    equ_ptrs[1]->assemble_bilinear_forms (local_cells, is_cell_at_bd);
  // multigroup iterations
  mg_iterations (local_cells, is_cell_at_bd, equ_ptrs);
}

template <int dim>
void MGBase<dim>::mg_iterations
(std::vector<typename DoFHandler<dim>::active_cell_iterator> &local_cells,
 std::vector<bool> &is_cell_at_bd,
 std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs)
{// this function needs to be overridden if JFNK is desired
  // by default, we give out Jacobi iteration scheme
  for (unsigned int g=0; g<n_group; <#increment#>) {
    <#statements#>
  }
  /*
  for (unsigned int g=0; g<n_group; ++g)
  {
    generate_group_rhses (sys_rhses, g);
    win_ptr->solve_in_group (sys_mats, g)
  }
   */
  // GS
  /*
  for (unsigned int g=0; g<n_group; ++g)
  {
    generate_group_rhses (sys_rhses, g);
    win_ptr->solve_in_group (sys_mats,vec_aflx,sys_rhses)
  }
   */
  // Jacobi
  /*
   for (unsigned int g=0; g<n_group; ++g)
     generate_group_rhses (sys_rhses, g);
   for (unsigned int g=0; g<n_group; ++g)
     win_ptr->solve_in_group (sys_mats,vec_aflx,sys_rhses)
   */
}

template class MGBase<2>;
template class MGBase<3>;
