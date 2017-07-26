#include <deal.II/fe/fe_values.h>

#include <boost/algorithm/string.hpp>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/cell_id.h>

#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/solver_bicgstab.h>

#include <algorithm>

#include "transport_base.h"
#include "bart_tools.h"

using namespace dealii;

template <int dim>
TransportBase<dim>::TransportBase
(ParameterHandler &prm,
 const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
 const std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr,
 const std_cxx11::shared_ptr<MaterialProperties> mat_ptr)
:
transport_model_name(prm.get("transport model")),
aq_name(prm.get("angular quadrature name")),
n_material(prm.get_integer("number of materials")),
discretization(prm.get("spatial discretization")),
n_group(prm.get_integer("number of groups")),
n_azi(prm.get_integer("angular quadrature order")),
is_eigen_problem(prm.get_bool("do eigenvalue calculations")),
do_nda(prm.get_bool("do NDA")),
do_print_sn_quad(prm.get_bool("do print angular quadrature info")),
have_reflective_bc(prm.get_bool("have reflective BC")),
p_order(prm.get_integer("finite element polynomial degree"))
{
  this->process_input (msh_ptr, aqd_ptr, mat_ptr);
  if (transport_model_name=="ep" && discretization=="dfem")
    c_penalty = 1.0 * p_order * (p_order + 1.0);
  sflx_proc.resize (n_group);
  sflx_proc_prev_gen.resize (n_group);
}

template <int dim>
TransportBase<dim>::~TransportBase ()
{
}

template <int dim>
void TransportBase<dim>::process_input (msh_ptr, aqd_ptr, mat_ptr)
{
  // mesh related
  {
    relative_position_to_id = msh_ptr->get_id_map ();
    if (have_reflective_bc)
      is_reflective_bc = msh_ptr->get_reflective_bc_map ();
  }
  
  // aq data related
  {
    n_total_ho_vars = aqd_ptr->get_n_total_ho_vars ();
    n_azi = aqd_ptr->get_sn_order ();
    n_dir = aqd_ptr->get_n_dir ();
    component_index = aqd_ptr->get_component_index_map ();
    inverse_component_index = aqd_ptr->get_inv_component_map ();
    wi = aqd_ptr->get_angular_weights ();
    omega_i = aqd_ptr->get_all_directions ();
    if (transport_model_name=="ep" && discretization=="dfem")
      tensor_norms = aqd_ptr->get_tensor_norms ();
    
    if (have_reflective_bc)
      reflective_direction_index = aqd_ptr->get_reflective_direction_index_map ();
  }
  
  // material related
  {
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
void TransportBase<dim>::initialize_assembly_related_objects
(DoFHandler<dim> &dof_handler,
 FE_Poly<TensorProductPolynomials<dim>,dim,dim>* fe)
{
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
  AssertThrow(do_nda==false,
              ExcMessage("Moments are generated only without NDA"));
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
void TransportBase<dim>::update_ho_moments_in_fiss
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
  double global_fiss_source = Utilities::MPI::sum (fiss_source, MPI_COMM_WORLD);
  return global_fiss_source;
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

// explicit instantiation to avoid linking error
template class TransportBase<2>;
template class TransportBase<3>;
