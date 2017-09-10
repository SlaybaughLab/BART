#include "eigen_base.h"
#include "mg_base.h"

template <int dim>
EigenBase<dim>::EigenBase (ParameterHandler &prm) : IterationBase<dim> (prm),
err_k_tol(1.0e-6),
err_phi_tol(1.0e-5),
keff(1.0)
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
  // assemble system matrices
  for (unsigned int i=0; i<equ_ptrs.size(); ++i)
    equ_ptrs[i]->assemble_bilinear_form ();
  
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
  for (unsigned int g=0; g<n_group; ++g)
    sflxes_proc[g] = 1.0;
  
  fission_source = equ_ptrs[0]->estimate_fiss_source (this->sflxes_proc);
}

// Override this function to do specific eigenvalue iteration as desired
template <int dim>
void EigenBase<dim>::eigen_iterations
(std::vector<Vector<double> > &sflxes_proc,
 std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs)
{
}

template <int dim>
void EigenBase<dim>::update_fiss_source_keff
(std::vector<Vector<double> > &sflxes_proc,
 std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs)
{
  keff_prev = keff;
  fiss_source_prev = fission_source;
  fiss_source = equ_ptrs[0]->estimate_fiss_source (sflxes_proc);
  keff = estimate_k (fiss_source, fiss_source_prev, keff_prev_gen);
}

template <int dim>
double EigenBase<dim>::estimate_k (double &fiss_source,
                                   double &fiss_source_prev,
                                   double &k_prev)
{
  return k_prev * fiss_source / fiss_source_prev;
}

template <int dim>
double EigenBase<dim>::estimate_k_diff (double &k, double &k_prev)
{
  return std::fabs (k - k_prev)/k;
}

template <int dim>
double EigenBase<dim>::get_keff ()
{
  return keff;
}

template class EigenBase<2>;
template class EigenBase<3>;
