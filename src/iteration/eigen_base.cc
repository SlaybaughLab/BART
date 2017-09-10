#include "eigen_base.h"
#include "mg_base.h"

template <int dim>
EigenBase<dim>::EigenBase (ParameterHandler &prm) : IterationBase<dim> (prm),
err_k_tol(1.0e-6),
err_phi_tol(1.0e-5)
{
  mg_ptr = build_mg_iterations (prm);
}

template <int dim>
EigenBase<dim>::~EigenBase ()
{
}

template <int dim>
EigenBase<dim>::do_iterations
void MGBase<dim>::do_iterations
(std::vector<Vector<double> > &sflxes_proc,
 std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs)
{
  if (this->do_nda)
    AssertThrow (equ_ptrs.size==2,
                 ExcMessage("There should be two equation pointers if do NDA"));
  // assemble system matrices for transport equation
  for (unsigned int i=0; i<equ_ptrs.size(); ++i)
    equ_ptrs[i]->assemble_bilinear_forms ();
  if (this->do_nda)
    equ_ptrs[1]->assemble_closure_bilinear_form (equ_ptrs[0]);
  // initialize fission process
  initialize_fiss_process (sflxes_proc, equ_ptrs);
  // perform eigenvalue iterations
  eigen_iterations (sflxes_proc, equ_ptrs);
}

template <int dim>
void EigenBase<dim>::initialize_fiss_process
(std::vector<Vector<double> > &sflxes_proc,
 std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs)
{
  // initialize sflxes with unit values
  equ_ptrs[0]->initialize_sflxes_proc (sflxes_proc);
  // calculate fission source based on initial scalar fluxes
  fiss_src = equ_ptrs[0]->estimate_fiss_src (sflxes_proc);
  // initialize keff
  keff = 1.0;
}

// Override this function to do specific eigenvalue iteration as desired
template <int dim>
void EigenBase<dim>::eigen_iterations
(std::vector<Vector<double> > &sflxes_proc,
 std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs)
{
}

template <int dim>
void EigenBase<dim>::update_prev_sflxes_fiss_src_keff
(std::vector<Vector<double> >&sflxes_proc)
{
  // update scalar fluxes from previous eigen iteration
  for (unsigned int g=0; g<this->n_group; ++g)
    sflxes_proc_prev_eigen[g] = sflxes_proc[g];
  // update fission source from previous eigen iteration
  fiss_src_prev = fiss_src;
  // update keff from previous eigen iteration
  keff_prev = keff;
}

template <int dim>
void EigenBase<dim>::calculate_fiss_src_keff
(std::vector<Vector<double> > &sflxes_proc,
 std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs)
{
  fiss_src = equ_ptrs[0]->estimate_fiss_src (sflxes_proc);
  keff = estimate_k ();
}

template <int dim>
double EigenBase<dim>::estimate_k ()
{
  return keff_prev * fiss_src / fiss_src_prev;
}

template <int dim>
double EigenBase<dim>::estimate_k_diff ()
{
  return std::fabs (keff - keff_prev) / keff;
}

template <int dim>
void EigenBase<dim>::get_keff (double &k)
{
  k = keff;
}

template class EigenBase<2>;
template class EigenBase<3>;
