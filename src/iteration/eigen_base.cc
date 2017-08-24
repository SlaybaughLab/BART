#include "eigen_base.h"

template <int dim>
EigenBase<dim>::EigenBase () : IterationBase<dim> (),
err_k_tol(1.0e-6),
err_phi_tol(1.0e-7),
err_phi_eigen_tol(1.0e-5),
keff(1.0)
{
}

template <int dim>
EigenBase<dim>::~EigenBase ()
{
}

template <int dim>
EigenBase<dim>::do_iterations
(ParameterHandler &prm,
 std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
 std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr,
 std_cxx11::shared_ptr<MaterialProperties> mat_ptr,
 std::vector<PETScWrappers::MPI::SparseMatrix*> sys_mats)
{
  this->initialize_equations (prm, msh_ptr, aqd_ptr, mat_ptr);
  eigen_iterations (msh_ptr, aqd_ptr, mat_ptr, sys_mats);
}

template <int dim>
void EigenBase<dim>::generate_system_matrices
(std::vector<PETScWrappers::MPI::SparseMatrix*> sys_mats)
{
  // EquationBase will be instantiated in InGroupBase
  mgs_ptr->generate_system_matrices (sys_mats);
}

template <int dim>
void EigenBase<dim>::initialize_fiss_process ()
{
  for (unsigned int g=0; g<n_group; ++g)
    this->sflx_proc[g] = 1.0;
  
  fission_source = this->trm_ptr->estimate_fiss_source (this->sflx_proc);
}

template <int dim>
double EigenBase<dim>::estimate_k (double &fiss_source,
                                   double &fiss_source_prev,
                                   double &k_prev)
{
  return k_prev * fiss_source / fiss_source_prev;
}

template <int dim>
double EigenBase<dim>::estimate_k_err (double &k, double &k_prev)
{
  return std::fabs (k - k_prev)/k;
}

template <int dim>
void EigenBase<dim>::eigen_iteration (double &keff,
                                     std_cxx11::shared_ptr<MGSolver<dim> > mgs_ptr)
{
}

template class EigenBase<2>;
template class EigenBase<3>;
