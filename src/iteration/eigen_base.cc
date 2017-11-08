#include "eigen_base.h"

template <int dim>
EigenBase<dim>::EigenBase (const ParameterHandler &prm)
:
IterationBase<dim>(prm),
err_k_tol(1.0e-6),
err_phi_tol(1.0e-5)
{
  sflxes_proc_prev_eigen.resize (this->n_group);
}

template <int dim>
EigenBase<dim>::~EigenBase ()
{
}

template <int dim>
void EigenBase<dim>::do_iterations
(std::vector<Vector<double> > &sflxes_proc,
 std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs,
 std_cxx11::shared_ptr<IGBase<dim> > ig_ptr,
 std_cxx11::shared_ptr<MGBase<dim> > mg_ptr)
{
  if (!this->do_nda)
  {
    // assemble system matrices for transport equation
    equ_ptrs.back()->assemble_bilinear_form ();
    // initialize fission process
    initialize_fiss_process (sflxes_proc, equ_ptrs.front());
    // perform eigenvalue iterations
    eigen_iterations (sflxes_proc, equ_ptrs, ig_ptr, mg_ptr);
  }
}

template <int dim>
void EigenBase<dim>::initialize_fiss_process
(std::vector<Vector<double> > &sflxes_proc,
 std_cxx11::shared_ptr<EquationBase<dim> > &equ_ptr)
{
  // calculate fission source based on initial scalar fluxes
  fiss_src = equ_ptr->estimate_fiss_src (sflxes_proc);
  // initialize keff
  keff = 1.0;
}

// Override this function to do specific eigenvalue iteration as desired
template <int dim>
void EigenBase<dim>::eigen_iterations
(std::vector<Vector<double> > &sflxes_proc,
 std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs,
 std_cxx11::shared_ptr<IGBase<dim> > ig_ptr,
 std_cxx11::shared_ptr<MGBase<dim> > mg_ptr)
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
 std_cxx11::shared_ptr<EquationBase<dim> > equ_ptr)
{
  fiss_src = equ_ptr->estimate_fiss_src (sflxes_proc);
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
double EigenBase<dim>::get_keff ()
{
  return keff;
}

template class EigenBase<2>;
template class EigenBase<3>;
