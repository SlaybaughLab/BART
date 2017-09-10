#include "eigen_base"

template <int dim>
PowerIteration<dim>::PowerIteration (ParameterHandler &)
:
EigenBase<dim> ()
{
}

template <int dim>
PowerIteration<dim>::~PowerIteration ()
{
}

template <int dim>
void PowerIteration<dim>::eigen_iterations
(std::vector<Vector<double> > &sflx_proc,
 std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs)
{
  double err_k = 1.0;
  double err_phi = 1.0;
  this->initialize_fiss_process (sflx_proc, equ_ptrs);
  while (err_k>this->err_k_tol || err_phi>this->err_phi_tol)
  {
    this->
    equ_ptrs[0]->scale_fiss_transfer_matrices (this->keff);
    equ_ptrs[0]->assemble_fixed_linear_form (sflx_proc);
    this->mg_ptr->mg_iterations (sflx_proc, equ_ptrs);
    this->update_fiss_source_keff ();
    err_phi = estimate_phi_diff (sflx_proc, this->sflx_proc_old);
    err_k = this->estimate_k_diff (this->keff, this->keff_prev);
    pout
    << "PI iter: " << ct << ", k: " << keff
    << ", err_k: " << err_k << ", err_phi: " << err_phi << std::endl;
  }
}

template class PowerIteration<2>;
template class PowerIteration<3>;
