template <int dim>
void PowerIteration<dim>::power_iteration
(std::vector<PETScWrappers::MPI::Vector*> &vec_ho_sflx)
{
  double err_k = 1.0;
  double err_phi = 1.0;
  unsigned int ct = 0;
  initialize_fiss_process (vec_ho_sflx);
  while (err_k>err_k_tol || err_phi>err_phi_eigen_tol)
  {
    ct += 1;
    update_ho_moments_in_fiss ();
    trm_ptr->scale_fiss_transfer_matrices ();//???
    trm_ptr->generate_ho_fixed_source (vec_ho_fixed_rhs,
                                       sflx_proc_prev_gen);
    source_iteration ();
    update_fiss_source_keff ();
    err_phi = estimate_phi_diff (vec_ho_sflx, vec_ho_sflx_prev_gen);
    err_k = std::fabs (keff - keff_prev_gen) / keff;
    pout
    << "PI iter: " << ct << ", k: " << keff
    << ", err_k: " << err_k << ", err_phi: " << err_phi << std::endl;
  }
}

template class PowerIteration<2>;
template class PowerIteration<3>;
