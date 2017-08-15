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
TransportBase<dim>::TransportBase (ParameterHandler &prm)
:
mpi_communicator (MPI_COMM_WORLD),
triangulation (mpi_communicator,
               typename Triangulation<dim>::MeshSmoothing
               (Triangulation<dim>::smoothing_on_refinement |
                Triangulation<dim>::smoothing_on_coarsening)),
dof_handler (triangulation),
err_k_tol(1.0e-6),
err_phi_tol(1.0e-7),
err_phi_eigen_tol(1.0e-5),
ho_linear_solver_name(prm.get("HO linear solver name")),
ho_preconditioner_name(prm.get("HO preconditioner name")),
pcout(std::cout,
      (Utilities::MPI::this_mpi_process(mpi_communicator)
       == 0))
{
  initialize_aq (prm);
  n_total_ho_vars = aqd_ptr->get_n_total_ho_vars ();
  sol_ptr = std_cxx11::shared_ptr<PreconditionerSolver>
  (new PreconditionerSolver (prm, n_total_ho_vars, mpi_communicator));
  def_ptr = std_cxx11::shared_ptr<ProblemDefinition>
  (new ProblemDefinition(prm));
  msh_ptr = std_cxx11::shared_ptr<MeshGenerator<dim> >
  (new MeshGenerator<dim>(prm));
  mat_ptr = std_cxx11::shared_ptr<MaterialProperties>
  (new MaterialProperties(prm));
  this->process_input ();
  sflx_proc.resize (n_group);
  sflx_proc_prev_gen.resize (n_group);
}

template <int dim>
TransportBase<dim>::~TransportBase ()
{
  dof_handler.clear();
}

template <int dim>
void TransportBase<dim>::process_input ()
{
  // basic parameters
  {
    // from basic problem definition
    transport_model_name = def_ptr->get_transport_model ();
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
void TransportBase<dim>::initialize_aq (ParameterHandler &prm)
{
  aq_name = prm.get ("angular quadrature name");
  AssertThrow (aq_name=="lsgc",
               ExcMessage("only LS-GC quadrature is implemented now."));
  if (aq_name=="lsgc")
    aqd_ptr = std_cxx11::shared_ptr<AQBase<dim> > (new AQLSGC<dim>(prm));
  aqd_ptr->make_aq (prm);
}

template <int dim>
void TransportBase<dim>::report_system ()
{
  pcout << "SN quadrature order: " << n_azi << std::endl
  << "Number of angles: " << n_dir << std::endl
  << "Number of groups: " << n_group << std::endl;

  radio ("Transport model", transport_model_name);
  radio ("Spatial discretization", discretization);
  radio ("HO linear solver", ho_linear_solver_name);
  if (ho_linear_solver_name!="direct")
    radio ("HO preconditioner", ho_preconditioner_name);
  radio ("do NDA?", do_nda);
  
  radio ("Number of cells", triangulation.n_global_active_cells());
  radio ("High-order total DoF counts", n_total_ho_vars*dof_handler.n_dofs());

  if (is_eigen_problem)
    radio ("Problem type: k-eigenvalue problem");
  if (do_nda)
    radio ("NDA total DoF counts", n_group*dof_handler.n_dofs());
  radio ("print sn quad?", do_print_sn_quad);
  if (do_print_sn_quad &&
      Utilities::MPI::this_mpi_process(mpi_communicator)==0)
    aqd_ptr->print_angular_quad ();
  radio ("is eigenvalue problem?", is_eigen_problem);
}

template <int dim>
void TransportBase<dim>::setup_system ()
{
  radio ("setup system");
  initialize_dealii_objects ();
  initialize_system_matrices_vectors ();
}

template <int dim>
void TransportBase<dim>::initialize_system_matrices_vectors ()
{
  DynamicSparsityPattern dsp (relevant_dofs);

  if (discretization=="dfem")
  {
    /*
     Table<2,DoFTools::Coupling> cell_coupling (1,1);
     Table<2,DoFTools::Coupling> face_coupling (1,1);

     cell_coupling[0][0] = DoFTools::nonzero;
     face_coupling[0][0] = DoFTools::nonzero;

     DoFTools::make_flux_sparsity_pattern (dof_handler,
     dsp,
     cell_coupling,
     face_coupling);
     */

    DoFTools::make_flux_sparsity_pattern (dof_handler,
                                          dsp,
                                          constraints,
                                          false);
  }
  else
    DoFTools::make_sparsity_pattern (dof_handler,
                                     dsp,
                                     constraints,
                                     false);

  // be careful with the following
  SparsityTools::distribute_sparsity_pattern (dsp,
                                              dof_handler.n_locally_owned_dofs_per_processor (),
                                              mpi_communicator,
                                              relevant_dofs);

  for (unsigned int g=0; g<n_group; ++g)
  {
    if (do_nda)
    {
      vec_lo_sys.push_back (new LA::MPI::SparseMatrix);
      vec_lo_rhs.push_back (new LA::MPI::Vector);
      vec_lo_sflx.push_back (new LA::MPI::Vector);
      vec_lo_sflx_old.push_back (new LA::MPI::Vector);
      vec_lo_fixed_rhs.push_back (new LA::MPI::Vector);
    }

    vec_ho_sflx.push_back (new LA::MPI::Vector);
    vec_ho_sflx_prev_gen.push_back (new LA::MPI::Vector);
    vec_ho_sflx_old.push_back (new LA::MPI::Vector);

    for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
    {
      vec_ho_sys.push_back (new LA::MPI::SparseMatrix);
      vec_aflx.push_back (new LA::MPI::Vector);
      vec_ho_rhs.push_back (new LA::MPI::Vector);
      vec_ho_fixed_rhs.push_back (new LA::MPI::Vector);
    }
  }

  for (unsigned int g=0; g<n_group; ++g)
  {
    if (do_nda)
    {
      vec_lo_sys[g]->reinit (local_dofs,
                             local_dofs,
                             dsp,
                             mpi_communicator);
      vec_lo_rhs[g]->reinit (local_dofs,
                             mpi_communicator);
      vec_lo_fixed_rhs[g]->reinit (local_dofs,
                                   mpi_communicator);
      vec_lo_sflx[g]->reinit (local_dofs,
                              mpi_communicator);
      vec_lo_sflx_old[g]->reinit (local_dofs,
                                  mpi_communicator);
    }

    vec_ho_sflx[g]->reinit (local_dofs,
                            mpi_communicator);
    vec_ho_sflx_old[g]->reinit (local_dofs,
                                mpi_communicator);

    for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
    {
      vec_ho_sys[get_component_index(i_dir, g)]->reinit(local_dofs,
                                                        local_dofs,
                                                        dsp,
                                                        mpi_communicator);
      vec_aflx[get_component_index(i_dir, g)]->reinit(local_dofs,
                                                      mpi_communicator);
      vec_ho_rhs[get_component_index(i_dir, g)]->reinit (local_dofs,
                                                         mpi_communicator);
      vec_ho_fixed_rhs[get_component_index(i_dir, g)]->reinit (local_dofs,
                                                               mpi_communicator);
    }
  }
}

template <int dim>
void TransportBase<dim>::assemble_ho_system ()
{
  radio ("Assemble volumetric bilinear forms");
  assemble_ho_volume_boundary ();

  if (discretization=="dfem")
  {
    AssertThrow (transport_model_name=="ep",
                 ExcMessage("DFEM is only implemented for even parity"));
    radio ("Assemble cell interface bilinear forms for DFEM");
    assemble_ho_interface ();
  }
}

template <int dim>
void TransportBase<dim>::initialize_dealii_objects ()
{
  if (discretization=="dfem")
    fe = (new FE_DGQ<dim> (p_order));
  else
    fe = (new FE_Q<dim> (p_order));

  dof_handler.distribute_dofs (*fe);

  local_dofs = dof_handler.locally_owned_dofs ();
  DoFTools::extract_locally_relevant_dofs (dof_handler,
                                           relevant_dofs);

  constraints.clear ();
  constraints.reinit (relevant_dofs);
  DoFTools::make_hanging_node_constraints (dof_handler,
                                           constraints);
  constraints.close ();

  q_rule = std_cxx11::shared_ptr<QGauss<dim> > (new QGauss<dim> (p_order + 1));
  qf_rule = std_cxx11::shared_ptr<QGauss<dim-1> > (new QGauss<dim-1> (p_order + 1));

  fv = std_cxx11::shared_ptr<FEValues<dim> >
  (new FEValues<dim> (*fe, *q_rule,
                      update_values | update_gradients |
                      update_quadrature_points |
                      update_JxW_values));

  fvf = std_cxx11::shared_ptr<FEFaceValues<dim> >
  (new FEFaceValues<dim> (*fe, *qf_rule,
                          update_values | update_gradients |
                          update_quadrature_points | update_normal_vectors |
                          update_JxW_values));

  fvf_nei = std_cxx11::shared_ptr<FEFaceValues<dim> >
  (new FEFaceValues<dim> (*fe, *qf_rule,
                          update_values | update_gradients |
                          update_quadrature_points | update_normal_vectors |
                          update_JxW_values));

  dofs_per_cell = fe->dofs_per_cell;
  n_q = q_rule->size();
  n_qf = qf_rule->size();

  local_dof_indices.resize (dofs_per_cell);
  neigh_dof_indices.resize (dofs_per_cell);
}

template <int dim>
void TransportBase<dim>::assemble_ho_volume_boundary ()
{
  // volumetric pre-assembly matrices
  std::vector<std::vector<FullMatrix<double> > >
  streaming_at_qp (n_q, std::vector<FullMatrix<double> > (n_dir, FullMatrix<double> (dofs_per_cell, dofs_per_cell)));

  std::vector<FullMatrix<double> >
  collision_at_qp (n_q, FullMatrix<double>(dofs_per_cell, dofs_per_cell));
  
  for (unsigned int i=0; i<local_cells.size(); ++i)
    vec_test_at_qp.push_back (FullMatrix<double> (n_q, dofs_per_cell));
  
  // this sector is for pre-assembling streaming and collision matrices at quadrature
  // points
  {
    typename DoFHandler<dim>::active_cell_iterator cell = local_cells[0];
    fv->reinit (cell);
    pre_assemble_cell_matrices (fv, cell, streaming_at_qp, collision_at_qp);
  }

  for (unsigned int k=0; k<n_total_ho_vars; ++k)
  {
    unsigned int g = get_component_group (k);
    unsigned int i_dir = get_component_direction (k);
    radio ("Assembling Component",k,"direction",i_dir,"group",g);
    FullMatrix<double> local_mat (dofs_per_cell, dofs_per_cell);

    for (unsigned int ic=0; ic<local_cells.size(); ++ic)
    {
      typename DoFHandler<dim>::active_cell_iterator cell = local_cells[ic];
      fv->reinit (cell);
      cell->get_dof_indices (local_dof_indices);
      local_mat = 0;
      integrate_cell_bilinear_form (fv,
                                    cell,
                                    local_mat,
                                    i_dir,
                                    g,
                                    streaming_at_qp,
                                    collision_at_qp);

      if (is_cell_at_bd[ic])
        for (unsigned int fn=0; fn<GeometryInfo<dim>::faces_per_cell; ++fn)
          if (cell->at_boundary(fn))
          {
            fvf->reinit (cell, fn);
            integrate_boundary_bilinear_form (fvf,
                                              cell,
                                              fn,
                                              local_mat,
                                              i_dir,
                                              g);
          }
      
      if (k==0)
        for (unsigned int qi=0; qi<n_q; ++qi)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            vec_test_at_qp[ic](qi, i) = fv->shape_value (i,qi) * fv->JxW (qi);
      
      vec_ho_sys[k]->add (local_dof_indices,
                          local_dof_indices,
                          local_mat);
    }
    vec_ho_sys[k]->compress (VectorOperation::add);
    pcout << "sys norm: " << vec_ho_sys[k]->l1_norm () << std::endl;
  }// components
}

// The following is a virtual function for integraing cell bilinear form;
// It can be overriden if cell pre-assembly is desirable
template <int dim>
void TransportBase<dim>::
pre_assemble_cell_matrices
(const std_cxx11::shared_ptr<FEValues<dim> > fv,
 typename DoFHandler<dim>::active_cell_iterator &cell,
 std::vector<std::vector<FullMatrix<double> > > &streaming_at_qp,
 std::vector<FullMatrix<double> > &collision_at_qp)
{// this is a virtual function
}

// The following is a virtual function for integraing cell bilinear form;
// It must be overriden
template <int dim>
void TransportBase<dim>::integrate_cell_bilinear_form
(const std_cxx11::shared_ptr<FEValues<dim> > fv,
 typename DoFHandler<dim>::active_cell_iterator &cell,
 FullMatrix<double> &cell_matrix,
 unsigned int &i_dir,
 unsigned int &g,
 std::vector<std::vector<FullMatrix<double> > > &streaming_at_qp,
 std::vector<FullMatrix<double> > &collision_at_qp)
{
}

// The following is a virtual function for integraing boundary bilinear form;
// It must be overriden
template <int dim>
void TransportBase<dim>::integrate_boundary_bilinear_form
(const std_cxx11::shared_ptr<FEFaceValues<dim> > fvf,
 typename DoFHandler<dim>::active_cell_iterator &cell,
 unsigned int &fn,/*face number*/
 FullMatrix<double> &cell_matrix,
 unsigned int &i_dir,
 unsigned int &g)
{// this is a virtual function
}

template <int dim>
void TransportBase<dim>::integrate_reflective_boundary_linear_form
(const std_cxx11::shared_ptr<FEFaceValues<dim> > fvf,
 typename DoFHandler<dim>::active_cell_iterator &cell,
 unsigned int &fn,/*face number*/
 std::vector<Vector<double> > &cell_rhses,
 unsigned int &i_dir,
 unsigned int &g)
{// this is a virtual function
}

template <int dim>
void TransportBase<dim>::assemble_ho_interface ()
{
  FullMatrix<double> vp_up (dofs_per_cell, dofs_per_cell);
  FullMatrix<double> vp_un (dofs_per_cell, dofs_per_cell);
  FullMatrix<double> vn_up (dofs_per_cell, dofs_per_cell);
  FullMatrix<double> vn_un (dofs_per_cell, dofs_per_cell);

  for (unsigned int k=0; k<n_total_ho_vars; ++k)
  {
    unsigned int g = get_component_group (k);
    unsigned int i_dir = get_component_direction (k);

    for (unsigned int ic=0; ic<local_cells.size(); ++ic)
    {
      typename DoFHandler<dim>::active_cell_iterator
      cell = local_cells[ic];
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int fn=0; fn<GeometryInfo<dim>::faces_per_cell; ++fn)
        if (!cell->at_boundary(fn) &&
            cell->neighbor(fn)->id()<cell->id())
        {
          fvf->reinit (cell, fn);
          typename DoFHandler<dim>::cell_iterator
          neigh = cell->neighbor(fn);
          neigh->get_dof_indices (neigh_dof_indices);
          fvf_nei->reinit (neigh, cell->neighbor_face_no(fn));

          vp_up = 0;
          vp_un = 0;
          vn_up = 0;
          vn_un = 0;

          integrate_interface_bilinear_form (fvf, fvf_nei,/*FEFaceValues objects*/
                                             cell, neigh,/*cell iterators*/
                                             fn,
                                             i_dir, g,/*specific component*/
                                             vp_up, vp_un, vn_up, vn_un);
          vec_ho_sys[k]->add (local_dof_indices,
                              local_dof_indices,
                              vp_up);

          vec_ho_sys[k]->add (local_dof_indices,
                              neigh_dof_indices,
                              vp_un);

          vec_ho_sys[k]->add (neigh_dof_indices,
                              local_dof_indices,
                              vn_up);

          vec_ho_sys[k]->add (neigh_dof_indices,
                              neigh_dof_indices,
                              vn_un);
        }// target faces
    }
    vec_ho_sys[k]->compress(VectorOperation::add);
  }// component
}

// The following is a virtual function for integrating DG interface for HO system
// it must be overriden
template <int dim>
void TransportBase<dim>::integrate_interface_bilinear_form
(const std_cxx11::shared_ptr<FEFaceValues<dim> > fvf,
 const std_cxx11::shared_ptr<FEFaceValues<dim> > fvf_nei,
 typename DoFHandler<dim>::active_cell_iterator &cell,
 typename DoFHandler<dim>::cell_iterator &neigh,/*cell iterator for cell*/
 unsigned int &fn,/*concerning face number in local cell*/
 unsigned int &i_dir,
 unsigned int &g,
 FullMatrix<double> &vp_up,
 FullMatrix<double> &vp_un,
 FullMatrix<double> &vn_up,
 FullMatrix<double> &vn_un)
{
}

template <int dim>
void TransportBase<dim>::generate_moments ()
{
  // FitIt: only scalar flux is generated for now
  AssertThrow(do_nda==false, ExcMessage("Moments are generated only without NDA"));
  if (!do_nda)
    for (unsigned int g=0; g<n_group; ++g)
    {
      *vec_ho_sflx_old[g] = *vec_ho_sflx[g];
      *vec_ho_sflx[g] = 0;
      for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
        vec_ho_sflx[g]->add (wi[i_dir], *vec_aflx[get_component_index(i_dir, g)]);
      sflx_proc[g] = *vec_ho_sflx[g];
    }
}

template <int dim>
void TransportBase<dim>::generate_ho_rhs ()
{
}

template <int dim>
void TransportBase<dim>::NDA_PI ()
{
}

template <int dim>
void TransportBase<dim>::NDA_SI ()
{
}

template <int dim>
void TransportBase<dim>::scale_fiss_transfer_matrices ()
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
void TransportBase<dim>::generate_ho_fixed_source ()
{
}

template <int dim>
void TransportBase<dim>::initialize_fiss_process ()
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
void TransportBase<dim>::update_ho_moments_in_fiss ()
{
  for (unsigned int g=0; g<n_group; ++g)
  {
    *vec_ho_sflx_prev_gen[g] = *vec_ho_sflx[g];
    sflx_proc_prev_gen[g] = *vec_ho_sflx_prev_gen[g];
  }
}

template <int dim>
void TransportBase<dim>::update_fiss_source_keff ()
{
  keff_prev_gen = keff;
  fission_source_prev_gen = fission_source;
  fission_source = estimate_fiss_source (sflx_proc);
  keff = estimate_k (fission_source, fission_source_prev_gen, keff_prev_gen);
  //renormalize_sflx (vec_ho_sflx);
}

template <int dim>
void TransportBase<dim>::power_iteration ()
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
void TransportBase<dim>::source_iteration ()
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
void TransportBase<dim>::renormalize_sflx
(std::vector<LA::MPI::Vector*> &target_sflxes)
{
  AssertThrow (target_sflxes.size()==n_group,
               ExcMessage("vector of scalar fluxes must have a size of n_group"));
  double norm_factor = target_sflxes[0]->max ();
  for (unsigned int g=0; g<n_group; ++g)
    *target_sflxes[g] /= norm_factor;
}

template <int dim>
void TransportBase<dim>::postprocess ()
{// do nothing in the base class
}

template <int dim>
double TransportBase<dim>::estimate_fiss_source (std::vector<Vector<double> > &phis_this_process)
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
double TransportBase<dim>::estimate_k (double &fiss_source,
                                       double &fiss_source_prev_gen,
                                       double &k_prev_gen)
{
  return k_prev_gen * fiss_source / fiss_source_prev_gen;
}

template <int dim>
double TransportBase<dim>::estimate_phi_diff
(std::vector<LA::MPI::Vector*> &phis_newer,
 std::vector<LA::MPI::Vector*> &phis_older)
{
  AssertThrow (phis_newer.size ()== phis_older.size (),
               ExcMessage ("n_groups for different phis should be identical"));
  double err = 0.0;
  for (unsigned int i=0; i<phis_newer.size (); ++i)
  {
    LA::MPI::Vector dif = *(phis_newer)[i];
    dif -= *(phis_older)[i];
    err = std::max (err, dif.l1_norm () / phis_newer[i]->l1_norm ());
  }
  return err;
}

template <int dim>
void TransportBase<dim>::do_iterations ()
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

template <int dim>
void TransportBase<dim>::output_results () const
{
  std::string sec_name = "Graphical output";
  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);

  for (unsigned int g=0; g<n_group; ++g)
  {
    std::ostringstream os;
    os << "ho_phi_g_" << g;
    data_out.add_data_vector (sflx_proc[g], os.str ());
  }

  Vector<float> subdomain (triangulation.n_active_cells ());
  for (unsigned int i=0; i<subdomain.size(); ++i)
    subdomain(i) = triangulation.locally_owned_subdomain ();
  data_out.add_data_vector (subdomain, "subdomain");

  data_out.build_patches ();

  const std::string filename =
  (namebase + "-" + discretization + "-" + Utilities::int_to_string
   (triangulation.locally_owned_subdomain (), 4));
  std::ofstream output ((filename + ".vtu").c_str ());
  data_out.write_vtu (output);

  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
  {
    std::vector<std::string> filenames;
    for (unsigned int i=0;
         i<Utilities::MPI::n_mpi_processes(mpi_communicator);
         ++i)
      filenames.push_back (namebase + "-" + discretization + "-" +
                           Utilities::int_to_string (i, 4) + ".vtu");
    std::ostringstream os;
    os << namebase << "-" << discretization << "-" << global_refinements << ".pvtu";
    std::ofstream master_output ((os.str()).c_str ());
    data_out.write_pvtu_record (master_output, filenames);
  }
}

template <int dim>
void TransportBase<dim>::run ()
{
  radio ("making grid");
  msh_ptr->make_grid (triangulation);
  msh_ptr->get_relevant_cell_iterators (dof_handler,
                                        local_cells,
                                        ref_bd_cells,
                                        is_cell_at_bd,
                                        is_cell_at_ref_bd);
  //msh_ptr.reset ();
  setup_system ();
  report_system ();
  assemble_ho_system ();
  do_iterations ();
  output_results();
}

// wrapper functions used to retrieve info from various Hash tables
template <int dim>
unsigned int TransportBase<dim>::get_component_index
(unsigned int incident_angle_index, unsigned int g)
{
  // retrieve component indecis given direction and group
  // must be used after initializing the index map
  return component_index[std::make_pair (incident_angle_index, g)];
}

template <int dim>
unsigned int TransportBase<dim>::get_component_direction (unsigned int comp_ind)
{
  return inverse_component_index[comp_ind].first;
}

template <int dim>
unsigned int TransportBase<dim>::get_component_group (unsigned int comp_ind)
{
  return inverse_component_index[comp_ind].second;
}

template <int dim>
unsigned int TransportBase<dim>::get_reflective_direction_index
(unsigned int boundary_id, unsigned int incident_angle_index)
{
  AssertThrow (is_reflective_bc[boundary_id],
               ExcMessage ("must be reflective boundary to retrieve the reflective boundary"));
  return reflective_direction_index[std::make_pair (boundary_id,
                                                    incident_angle_index)];
}

//functions used to cout information for diagonose or just simply cout
template <int dim>
void TransportBase<dim>::radio (std::string str)
{
  pcout << str << std::endl;
}

template <int dim>
void TransportBase<dim>::radio (std::string str1, std::string str2)
{
  pcout << str1 << ": " << str2 << std::endl;
}

template <int dim>
void TransportBase<dim>::radio (std::string str,
                                double num)
{
  pcout << str << ": " << num << std::endl;
}

template <int dim>
void TransportBase<dim>::radio (std::string str1, unsigned int num1,
                                std::string str2, unsigned int num2,
                                std::string str3, unsigned int num3)
{
  pcout << str1 << ": " << num1 << ", ";
  pcout << str2 << ": " << num2 << ", ";
  pcout << str3 << ": " << num3 << std::endl;;
}

template <int dim>
void TransportBase<dim>::radio (std::string str,
                                unsigned int num)
{
  pcout << str << ": " << num << std::endl;
}

template <int dim>
void TransportBase<dim>::radio (std::string str, bool boolean)
{
  pcout << str << ": " << (boolean?"true":"false") << std::endl;
}

template <int dim>
void TransportBase<dim>::radio ()
{
  pcout << "-------------------------------------" << std::endl << std::endl;
}

// explicit instantiation to avoid linking error
template class TransportBase<2>;
template class TransportBase<3>;
