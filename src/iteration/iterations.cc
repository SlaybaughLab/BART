#include <deal.II/base/index_set.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/cell_id.h>

#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/solver_bicgstab.h>

#include <algorithm>

#include "iterations.h"
#include "../aqdata/aq_base.h"
#include "../aqdata/aq_lsgc.h"

using namespace dealii;

template <int dim>
Iterations<dim>::Iterations
(ParameterHandler &prm,
 const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
 const std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr)
:
is_eigen_problem(prm.get_bool("do eigenvalue calculations")),
do_nda(prm.get_bool("do NDA")),
n_group(prm.get_integer("number of groups"))
{
  n_total_ho_vars = aqd_ptr->get_n_total_ho_vars ();
  sol_ptr = build_solution (prm, n_total_ho_vars);
  mat_ptr = build_material (prm);
  msh_ptr->get_relevant_cell_iterators (dof_handler,
                                        local_cells,
                                        ref_bd_cells,
                                        is_cell_at_bd,
                                        is_cell_at_ref_bd);
  trm_ptr = build_transport_model ();
  sflx_proc.resize (n_group);
  sflx_proc_prev_gen.resize (n_group);
}

template <int dim>
Iterations<dim>::~Iterations ()
{
}

template <int dim>
void Iterations<dim>::initialize_system_matrices_vectors
(SparsityPatternType &dsp, IndexSet &local_dofs)
{
  for (unsigned int g=0; g<n_group; ++g)
  {
    if (do_nda)
    {
      vec_lo_sys.push_back (new LA::MPI::SparseMatrix);
      //vec_lo_rhs.push_back (new LA::MPI::Vector);
      vec_lo_sflx.push_back (new LA::MPI::Vector);
      //vec_lo_sflx_old.push_back (new LA::MPI::Vector);
      //vec_lo_fixed_rhs.push_back (new LA::MPI::Vector);
    }
    
    vec_ho_sflx.push_back (new LA::MPI::Vector);
    //vec_ho_sflx_prev_gen.push_back (new LA::MPI::Vector);
    //vec_ho_sflx_old.push_back (new LA::MPI::Vector);
    
    for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
    {
      vec_ho_sys.push_back (new LA::MPI::SparseMatrix);
      //vec_aflx.push_back (new LA::MPI::Vector);
      //vec_ho_rhs.push_back (new LA::MPI::Vector);
      //vec_ho_fixed_rhs.push_back (new LA::MPI::Vector);
    }
  }
  
  for (unsigned int g=0; g<n_group; ++g)
  {
    if (do_nda)
    {
      vec_lo_sys[g]->reinit (local_dofs,
                             local_dofs,
                             dsp,
                             MPI_COMM_WORLD);
      //vec_lo_rhs[g]->reinit (local_dofs, MPI_COMM_WORLD);
      //vec_lo_fixed_rhs[g]->reinit (local_dofs, MPI_COMM_WORLD);
      vec_lo_sflx[g]->reinit (local_dofs, MPI_COMM_WORLD);
      //vec_lo_sflx_old[g]->reinit (local_dofs, MPI_COMM_WORLD);
    }
    
    vec_ho_sflx[g]->reinit (local_dofs, MPI_COMM_WORLD);
    //vec_ho_sflx_old[g]->reinit (local_dofs, MPI_COMM_WORLD);
  }
  
  for (unsigned int k=0; k<n_total_ho_vars; ++k)
  {
    vec_ho_sys[k]->reinit(local_dofs,
                          local_dofs,
                          dsp,
                          MPI_COMM_WORLD);
    //vec_aflx[k]->reinit(local_dofs, MPI_COMM_WORLD);
    //vec_ho_rhs[k]->reinit (local_dofs, MPI_COMM_WORLD);
    //vec_ho_fixed_rhs[k]->reinit (local_dofs, MPI_COMM_WORLD);
  }
}

template <int dim>
void Iterations<dim>::scale_fiss_transfer_matrices ()
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
void Iterations<dim>::

template <int dim>
void Iterations<dim>::update_ho_moments_in_fiss
(std::vector<PETScWrappers::MPI::Vector*> &vec_ho_sflx,
 std::vector<PETScWrappers::MPI::Vector*> &vec_ho_sflx_prev_gen)
{
  for (unsigned int g=0; g<n_group; ++g)
  {
    *vec_ho_sflx_prev_gen[g] = *vec_ho_sflx[g];
    sflx_proc_prev_gen[g] = *vec_ho_sflx_prev_gen[g];
  }
}

template <int dim>
void Iterations<dim>::update_fiss_source_keff ()
{
  keff_prev_gen = keff;
  fission_source_prev_gen = fission_source;
  fission_source = trm_ptr->estimate_fiss_source (sflx_proc);
  keff = estimate_k (fission_source, fission_source_prev_gen, keff_prev_gen);
}

template <int dim>
void Iterations<dim>::power_iteration
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
    radio ();
  }
}

template <int dim>
void Iterations<dim>::solve_problems
(std::vector<Vector<double> > &sflx_this_proc)
{
  if (is_eigen_problem)
  {
    std_cxx11::shared_ptr<EigenBase<dim> > pro_ptr = build_eigen_problem (prm);
    pro_ptr->do_iterations ();
    pro_ptr->get_sflx_proc (sflx_this_proc);
    keff = pro_ptr->get_keff ();
  }
  else
  {
    std_cxx11::shared_ptr<MGBase<dim> > pro_ptr = build_mg_problem (prm);
    pro_ptr->do_iterations ();
    pro_ptr->get_sflx_proc (sflx_this_proc);
  }
}

template <int dim>
void Iterations<dim>::get_keff (double &keff)
{
  keff = this->keff;
}

// explicit instantiation to avoid linking error
template class Iterations<2>;
template class Iterations<3>;
