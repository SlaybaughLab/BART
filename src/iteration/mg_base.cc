#include "mg_base.h"

template <int dim>
MGBase<dim>::MGBase (ParameterHandler &prm)
:
IterationBase<dim> (prm),
err_phi_tol(1.0e-5)
{
  // TODO: needs change if anisotropic scattering is needed
  sflxes_proc_prev_mg.resize (this->n_group);
  ig_ptr = build_ig_iteration (prm);
}

template <int dim>
MGBase<dim>::~MGBase ()
{
}

template <int dim>
void MGBase<dim>::do_iterations
(std::vector<Vector<double> > &sflxes_proc,
 std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs)
{
  if (this->do_nda)
    AssertThrow (equ_ptrs.size==2,
                 ExcMessage("There should be two equation pointers if do NDA"));
  // assemble bilinear forms of available equations
  for (unsigned int i=0; i<equ_ptrs.size(); ++i)
    equ_ptrs[i]->assemble_bilinear_forms ();
  if (this->do_nda)
    equ_ptrs[1]->assemble_closure_bilinear_form (equ_ptrs[0]);
  // multigroup iterations
  mg_iterations (sflx_proc, equ_ptrs);
}

// virtual function for all multigroup iteration method. It has to be overriden
// per derived class of MGBase. If it's fixed source problem, it will be called
// internally in do_iterations. Otherwise, it will be called externally in EigenBase
// instances.
template <int dim>
void MGBase<dim>::mg_iterations
(std::vector<Vector<double> > &sflx_proc,
 std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs)
{// this function needs to be overridden if JFNK is desired
  // by default, we give out Jacobi iteration scheme
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
