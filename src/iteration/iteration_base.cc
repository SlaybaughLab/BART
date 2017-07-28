#include <deal.II/fe/fe_values.h>

#include <boost/algorithm/string.hpp>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/cell_id.h>

#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/solver_bicgstab.h>

#include <algorithm>

#include "transport_base.h"
#include "../aqdata/aq_base.h"
#include "../aqdata/aq_lsgc.h"

using namespace dealii;

template <int dim>
IterationBase<dim>::IterationBase
(ParameterHandler &prm,
 const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
 const std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr)
:
err_k_tol(1.0e-6),
err_phi_tol(1.0e-7),
err_phi_eigen_tol(1.0e-5),
ho_linear_solver_name(prm.get("HO linear solver name")),
ho_preconditioner_name(prm.get("HO preconditioner name")),
is_eigen_problem(prm.get_bool("do eigenvalue calculations")),
do_nda(prm.get_bool("do NDA")),
pcout(std::cout,
      Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
{
  n_total_ho_vars = aqd_ptr->get_n_total_ho_vars ();
  sol_ptr = build_solution (prm, n_total_ho_vars);
  mat_ptr = build_material (prm);
  msh_ptr->get_relevant_cell_iterators (dof_handler,
                                        local_cells,
                                        ref_bd_cells,
                                        is_cell_at_bd,
                                        is_cell_at_ref_bd);
  this->process_input ();
  sflx_proc.resize (n_group);
  sflx_proc_prev_gen.resize (n_group);
}

template <int dim>
IterationBase<dim>::~IterationBase ()
{
}

template <int dim>
void IterationBase<dim>::process_input ()
{
  // basic parameters
  {
    // from basic problem definition
    n_group = def_ptr->get_n_group ();
    n_material = mat_ptr->get_n_material ();
    p_order = def_ptr->get_fe_order ();
    discretization = def_ptr->get_discretization ();
    have_reflective_bc = def_ptr->get_reflective_bool ();
    do_nda = def_ptr->get_nda_bool ();
    is_eigen_problem = def_ptr->get_eigen_problem_bool ();
    do_print_sn_quad = def_ptr->get_print_sn_quad_bool ();
    global_refinements = def_ptr->get_uniform_refinement ();
    namebase = def_ptr->get_output_namebase ();

    // from angular quadrature data
    n_azi = aqd_ptr->get_sn_order ();
    n_dir = aqd_ptr->get_n_dir ();
    component_index = aqd_ptr->get_component_index_map ();
    inverse_component_index = aqd_ptr->get_inv_component_map ();
    wi = aqd_ptr->get_angular_weights ();
    omega_i = aqd_ptr->get_all_directions ();
    if (transport_model_name=="ep" &&
        discretization=="dfem")
    {
      tensor_norms = aqd_ptr->get_tensor_norms ();
      c_penalty = 1.0 * p_order * (p_order + 1.0);
    }
  }

  if (have_reflective_bc)
  {
    is_reflective_bc = msh_ptr->get_reflective_bc_map ();
    reflective_direction_index = aqd_ptr->get_reflective_direction_index_map ();
  }

  // material properties
  {
    relative_position_to_id = msh_ptr->get_id_map ();
    all_sigt = mat_ptr->get_sigma_t ();
    all_inv_sigt = mat_ptr->get_inv_sigma_t ();
    all_sigs = mat_ptr->get_sigma_s ();
    all_sigs_per_ster = mat_ptr->get_sigma_s_per_ster ();
    if (is_eigen_problem)
    {
      is_material_fissile = mat_ptr->get_fissile_id_map ();
      all_nusigf = mat_ptr->get_nusigf ();
      all_ksi_nusigf = mat_ptr->get_ksi_nusigf ();
      all_ksi_nusigf_per_ster = mat_ptr->get_ksi_nusigf_per_ster ();
    }
    else
    {
      all_q = mat_ptr->get_q ();
      all_q_per_ster = mat_ptr->get_q_per_ster ();
    }
  }
}

template <int dim>
void IterationBase<dim>::NDA_PI ()
{
}

template <int dim>
void IterationBase<dim>::NDA_SI ()
{
}

template <int dim>
void IterationBase<dim>::scale_fiss_transfer_matrices ()
{
  AssertThrow (do_nda==false,
               ExcMessage("we don't scale fission transfer without NDA"));
  if (!do_nda)
  {
    scaled_fiss_transfer_per_ster.resize (n_material);
    for (unsigned int m=0; m<n_material; ++m)
    {
      std::vector<std::vector<double> >  tmp (n_group, std::vector<double>(n_group));
      if (is_material_fissile[m])
        for (unsigned int gin=0; gin<n_group; ++gin)
          for (unsigned int g=0; g<n_group; ++g)
            tmp[gin][g] = all_ksi_nusigf_per_ster[m][gin][g] / keff;
      scaled_fiss_transfer_per_ster[m] = tmp;
    }
  }
}

template <int dim>
void IterationBase<dim>::initialize_fiss_process ()
{
  for (unsigned int g=0; g<n_group; ++g)
  {
    *vec_ho_sflx[g] = 1.0;
    sflx_proc[g] = *vec_ho_sflx[g];
  }
  fission_source = estimate_fiss_source (sflx_proc);
  keff = 1.0;
}

template <int dim>
void IterationBase<dim>::update_ho_moments_in_fiss ()
{
  for (unsigned int g=0; g<n_group; ++g)
  {
    *vec_ho_sflx_prev_gen[g] = *vec_ho_sflx[g];
    sflx_proc_prev_gen[g] = *vec_ho_sflx_prev_gen[g];
  }
}

template <int dim>
void IterationBase<dim>::update_fiss_source_keff ()
{
  keff_prev_gen = keff;
  fission_source_prev_gen = fission_source;
  fission_source = estimate_fiss_source (sflx_proc);
  keff = estimate_k (fission_source, fission_source_prev_gen, keff_prev_gen);
  //renormalize_sflx (vec_ho_sflx);
}

template <int dim>
void IterationBase<dim>::power_iteration ()
{
  double err_k = 1.0;
  double err_phi = 1.0;
  unsigned int ct = 0;
  initialize_fiss_process ();
  while (err_k>err_k_tol || err_phi>err_phi_eigen_tol)
  {
    ct += 1;
    update_ho_moments_in_fiss ();
    scale_fiss_transfer_matrices ();
    generate_ho_fixed_source ();
    source_iteration ();
    update_fiss_source_keff ();
    err_phi = estimate_phi_diff (vec_ho_sflx, vec_ho_sflx_prev_gen);
    err_k = std::fabs (keff - keff_prev_gen) / keff;
    pcout
    << "PI iter: " << ct << ", k: " << keff
    << ", err_k: " << err_k << ", err_phi: " << err_phi << std::endl;
    radio ();
  }
}

template <int dim>
void IterationBase<dim>::source_iteration ()
{
  unsigned int ct = 0;
  double err_phi = 1.0;
  double err_phi_old;
  //generate_moments ();
  while (err_phi>err_phi_tol)
  {
    //generate_ho_source ();
    ct += 1;
    generate_ho_rhs ();
    sol_ptr->ho_solve (vec_ho_sys,
                       vec_aflx,
                       vec_ho_rhs);
    generate_moments ();
    err_phi_old = err_phi;
    err_phi = estimate_phi_diff (vec_ho_sflx, vec_ho_sflx_old);
    double spectral_radius = err_phi / err_phi_old;
    pcout
    << "SI iter: " << ct
    << ", phi err: " << err_phi
    << ", spec. rad.: " << spectral_radius << std::endl;
  }
}

template <int dim>
void IterationBase<dim>::postprocess ()
{// do nothing in the base class
}

template <int dim>
double IterationBase<dim>::estimate_fiss_source (std::vector<Vector<double> > &phis_this_process)
{
  double fiss_source = 0.0;
  for (unsigned int ic=0; ic<local_cells.size(); ++ic)
  {
    typename DoFHandler<dim>::active_cell_iterator cell = local_cells[ic];
    std::vector<std::vector<double> > local_phis (n_group,
                                                  std::vector<double> (n_q));
    unsigned int material_id = cell->material_id ();
    if (is_material_fissile[material_id])
    {
      fv->reinit (cell);
      for (unsigned int g=0; g<n_group; ++g)
        fv->get_function_values (phis_this_process[g],
                                 local_phis[g]);
      for (unsigned int qi=0; qi<n_q; ++qi)
        for (unsigned int g=0; g<n_group; ++g)
          fiss_source += (all_nusigf[material_id][g] *
                          local_phis[g][qi] *
                          fv->JxW(qi));
    }
  }
  double global_fiss_source = Utilities::MPI::sum (fiss_source, mpi_communicator);
  return global_fiss_source;
}

template <int dim>
double IterationBase<dim>::estimate_k (double &fiss_source,
                                       double &fiss_source_prev_gen,
                                       double &k_prev_gen)
{
  return k_prev_gen * fiss_source / fiss_source_prev_gen;
}

template <int dim>
double IterationBase<dim>::estimate_phi_diff
(std::vector<PETScWrappers::MPI::Vector*> &phis_newer,
 std::vector<PETScWrappers::MPI::Vector*> &phis_older)
{
  AssertThrow (phis_newer.size ()== phis_older.size (),
               ExcMessage ("n_groups for different phis should be identical"));
  double err = 0.0;
  for (unsigned int i=0; i<phis_newer.size (); ++i)
  {
    PETScWrappers::MPI::Vector dif = *(phis_newer)[i];
    dif -= *(phis_older)[i];
    err = std::max (err, dif.l1_norm () / phis_newer[i]->l1_norm ());
  }
  return err;
}

template <int dim>
void IterationBase<dim>::do_iterations ()
{
  sol_ptr->initialize_ho_preconditioners (vec_ho_sys, vec_ho_rhs);
  if (is_eigen_problem)
  {
    if (do_nda)
      NDA_PI ();
    else
    {
      power_iteration ();
      postprocess ();
    }
  }
  else
  {
    if (do_nda)
      NDA_SI ();
    else
    {
      generate_ho_fixed_source ();
      generate_moments ();
      source_iteration ();
      postprocess ();
    }
  }
}

// explicit instantiation to avoid linking error
template class IterationBase<2>;
template class IterationBase<3>;
