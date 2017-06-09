#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_values.h>

#include <boost/algorithm/string.hpp>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/cell_id.h>

#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/solver_bicgstab.h>

#include <algorithm>

#include "../include/EP_SN.h"

using namespace dealii;

template <int dim>
EP_SN<dim>::EP_SN (ParameterHandler &prm)
:
ProblemDefinition<dim>(prm),
mpi_communicator (MPI_COMM_WORLD),
triangulation (mpi_communicator,
               typename Triangulation<dim>::MeshSmoothing
               (Triangulation<dim>::smoothing_on_refinement |
                Triangulation<dim>::smoothing_on_coarsening)),
dof_handler (triangulation),
err_k_tol(1.0e-6),
err_phi_tol(1.0e-6),
pcout(std::cout,
      (Utilities::MPI::this_mpi_process(mpi_communicator)
       == 0))
{
  this->process_input ();
  sflx_this_processor.resize (n_group);
}

template <int dim>
EP_SN<dim>::~EP_SN()
{
  dof_handler.clear();
  delete fe;
}

template <int dim>
void EP_SN<dim>::process_input ()
{
  {
    n_azi = this->get_sn_order ();
    global_refinements = this->get_uniform_refinement ();
    n_group = this->get_n_group ();
    n_dir = this->get_n_dir ();
    p_order = this->get_fe_order ();
    discretization = this->get_discretization ();
    pcout << "method " << discretization << std::endl;
    axis_max_values = this->get_axis_maxes ();
    ncell_per_dir = this->get_ncells ();
    
    relative_position_to_id = this->get_id_map ();
    cell_size_all_dir = this->get_cell_sizes ();
    have_reflective_bc = this->get_reflective_bool ();
    do_nda = this->get_nda_bool ();
    is_eigen_problem = this->get_eigen_problem_bool ();
    n_total_ho_vars = this->get_n_total_ho_vars ();
    do_print_sn_quad = this->get_print_sn_quad_bool ();
    if (do_print_sn_quad &&
        Utilities::MPI::this_mpi_process(mpi_communicator)==0)
      this->print_angular_quad ();
    
    component_index = this->get_component_index_map ();
    inverse_component_index = this->get_inv_component_map ();
    wi = this->get_angular_weights ();
    omega_i = this->get_all_directions ();
    tensor_norms = this->get_tensor_norms ();
  }
  
  if (have_reflective_bc)
  {
    is_reflective_bc = this->get_reflective_bc_map ();
    reflective_direction_index = this->get_reflective_direction_index_map ();
  }
  
  {
    relative_position_to_id = this->get_id_map ();
    all_sigt = this->get_sigma_t ();
    all_inv_sigt = this->get_inv_sigma_t ();
    all_sigs = this->get_sigma_s ();
    all_sigs_per_ster = this->get_sigma_s_per_ster ();
    if (is_eigen_problem)
    {
      is_material_fissile = this->get_fissile_id_map ();
      all_nusigf = this->get_nusigf ();
      all_ksi_nusigf = this->get_ksi_nusigf ();
      all_ksi_nusigf_per_ster = this->get_ksi_nusigf_per_ster ();
    }
    else
    {
      all_q = this->get_q ();
      all_q_per_ster = this->get_q_per_ster ();
    }
  }
}

template <int dim>
void EP_SN<dim>::get_cell_relative_position (Point<dim> &center,
                                             std::vector<unsigned int> &relative_position)
{
  AssertThrow (relative_position.size()==3,
               ExcMessage("relative position should be size 3 for any dimension"));
  if (dim>=1)
  {
    relative_position[0] = static_cast<unsigned int>(center[0] / cell_size_all_dir[0]);
    if (dim>=2)
    {
      relative_position[1] = static_cast<unsigned int>(center[1] / cell_size_all_dir[1]);
      if (dim==3)
        relative_position[2] = static_cast<unsigned int>(center[2] / cell_size_all_dir[2]);
    }
  }
}

template <int dim>
void EP_SN<dim>::get_cell_mfps (unsigned int &material_id, double &cell_dimension,
                                std::vector<double> &local_mfps)
{
  // estimate mean free path for input cell aiming for penalty coefficients
  // FixIt: find a better way to estimate
  AssertThrow (local_mfps.size()==n_group,
               ExcMessage("size of mfp should be identical to n_group"));
  for (unsigned int g=0; g<n_group; ++g)
    local_mfps[g] = all_sigt[material_id][g] * cell_dimension;
}


template <int dim>
unsigned int EP_SN<dim>::get_component_index (unsigned int &incident_angle_index,
                                              unsigned int &g)
{
  // retrieve component indecis given direction and group
  // must be used after initializing the index map
  return component_index[std::make_pair (incident_angle_index, g)];
}

template <int dim>
unsigned int EP_SN<dim>::get_direction (unsigned int &comp_ind)
{
  return inverse_component_index[comp_ind].first;
}

template <int dim>
unsigned int EP_SN<dim>::get_component_group (unsigned int &comp_ind)
{
  return inverse_component_index[comp_ind].second;
}

template <int dim>
unsigned int EP_SN<dim>::get_reflective_direction_index (unsigned int &boundary_id,
                                                         unsigned int &incident_angle_index)
{
  AssertThrow (is_reflective_bc[boundary_id],
               ExcMessage ("must be reflective boundary to retrieve the reflective boundary"));
  return reflective_direction_index[std::make_pair (boundary_id, incident_angle_index)];
}

template <int dim>
void EP_SN<dim>::generate_globally_refined_grid ()
{
  pcout << "generate refined grid" << std::endl;
  Point<dim> origin;
  Point<dim> diagonal;
  switch (dim)
  {
    case 1:
    {
      diagonal[0] = axis_max_values[0];
      break;
    }
      
    case 2:
    {
      diagonal[0] = axis_max_values[0];
      diagonal[1] = axis_max_values[1];
      break;
    }
      
    case 3:
    {
      diagonal[0] = axis_max_values[0];
      diagonal[1] = axis_max_values[1];
      diagonal[2] = axis_max_values[2];
      break;
    }
      
    default:
      break;
  }
  GridGenerator::subdivided_hyper_rectangle (triangulation,
                                             ncell_per_dir,
                                             origin,
                                             diagonal);
  triangulation.refine_global (global_refinements);
  pcout << "generate refined grid finished" << std::endl;
}

template <int dim>
void EP_SN<dim>::initialize_material_id ()
{
  for (typename Triangulation<dim>::active_cell_iterator
       cell=triangulation.begin_active();
       cell!=triangulation.end();
       ++cell)
    if (cell->is_locally_owned())
    {
      Point<dim> center = cell->center ();
      std::vector<unsigned int> relative_position (3);
      get_cell_relative_position (center, relative_position);
      unsigned int material_id = relative_position_to_id[relative_position];
      cell->set_material_id (material_id);
    }
}

template <int dim>
void EP_SN<dim>::report_system ()
{
  pcout << "SN quadrature order: "
  << n_azi
  << std::endl
  << "Number of angles: "
  << n_dir
  << std::endl
  << "Number of groups: "
  << n_group
  << std::endl;
  
  pcout << "Number of active cells: "
  << triangulation.n_global_active_cells()
  << std::endl
  << "Number of high-order degrees of freedom: "
  << n_total_ho_vars * dof_handler.n_dofs()
  << std::endl;
  
  if (is_eigen_problem)
    pcout << "Problem type: k-eigenvalue problem" << std::endl;
  
  if (do_nda)
    pcout << "NDA DoFs: "
    << n_group * dof_handler.n_dofs() * n_group
    << std::endl;
}

template <int dim>
void EP_SN<dim>::setup_system ()
{
  local_radio ("setup system");
  //TimerOutput::Scope t(computing_timer, "setup HO system");
  
  if (boost::iequals(discretization,"DFEM") || boost::iequals(discretization,"DG"))
    fe = new FE_DGQ<dim> (p_order);
  else
    fe = new FE_Q<dim> (p_order);
  
  dof_handler.distribute_dofs (*fe);
  
  local_dofs = dof_handler.locally_owned_dofs ();
  DoFTools::extract_locally_relevant_dofs (dof_handler,
                                           relevant_dofs);
  
  constraints.clear ();
  constraints.reinit (relevant_dofs);
  DoFTools::make_hanging_node_constraints (dof_handler,
                                           constraints);
  constraints.close ();
  
  local_radio ("dsp setup");
  DynamicSparsityPattern dsp (relevant_dofs);
  
  if (boost::iequals(discretization,"DFEM") ||
      boost::iequals(discretization,"DG"))
    DoFTools::make_flux_sparsity_pattern (dof_handler,
                                          dsp,
                                          constraints,
                                          false);
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
  
  local_radio("initialize system mats and vecs");
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
    vec_ho_sflx_old.push_back (new LA::MPI::Vector);
    vec_ho_fixed_rhs.push_back (new LA::MPI::Vector);
    
    for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
    {
      vec_ho_sys.push_back (new LA::MPI::SparseMatrix);
      vec_aflx.push_back (new LA::MPI::Vector);
      vec_ho_rhs.push_back (new LA::MPI::Vector);
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
    
    vec_ho_fixed_rhs[g]->reinit (local_dofs,
                                 mpi_communicator);
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
    }
  }
  
  c_penalty = 1.0 * p_order * (p_order + 1.0);
}

template <int dim>
void EP_SN<dim>::setup_boundary_ids ()
{
  AssertThrow (axis_max_values.size()==dim,
               ExcMessage("number of entries axis max values should be dimension"));
  
  typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active (), endc = triangulation.end ();
  for (; cell!=endc; ++cell)
  {
    if (cell->is_locally_owned())
    {
      for (unsigned int fn=0; fn<GeometryInfo<dim>::faces_per_cell; ++fn)
      {
        if (cell->face(fn)->at_boundary())
        {
          Point<dim> ct = cell->face(fn)->center();
          // left boundary
          if (std::fabs(ct[0])<1.0e-14)
          {
            cell->face(fn)->set_boundary_id (0);
          }
          
          // right boundary
          if (std::fabs(ct[0]-axis_max_values[0])<1.0e-14)
          {
            cell->face(fn)->set_boundary_id (1);
          }
          
          // 2D and 3D boundaries
          if (dim>1)
          {
            // 2D boundaries
            // front boundary
            if (std::fabs(ct[1])<1.0e-14)
              cell->face(fn)->set_boundary_id (2);
            
            // rear boundary
            if (std::fabs(ct[1]-axis_max_values[1])<1.0e-14)
              cell->face(fn)->set_boundary_id (3);
            
            // 3D boundaries
            if (dim>2)
            {
              // front boundary
              if (std::fabs(ct[2])<1.0e-14)
                cell->face(fn)->set_boundary_id (4);
              
              // rear boundary
              if (std::fabs(ct[2]-axis_max_values[2])<1.0e-14)
                cell->face(fn)->set_boundary_id (5);
            }
          }
        }
      }// face
    }// locally owned cell
  }// cell
}

template <int dim>
void EP_SN<dim>::assemble_ho_system ()
{
  local_radio ("Assemble volumetric bilinear forms");
  assemble_ho_volume_boundary_new ();
  
  if (boost::iequals(discretization,"DFEM") || boost::iequals(discretization,"DG"))
  {
    local_radio ("Assemble cell interface bilinear forms for DFEM");
    assemble_ho_interface_new ();
  }
}

template <int dim>
void EP_SN<dim>::assemble_ho_volume_boundary_new ()
{
  const QGauss<dim>  q_rule(p_order+1);
  const QGauss<dim-1>  qf_rule(p_order+1);
  
  // cell finite element object
  FEValues<dim> fv(*fe, q_rule,
                   update_values | update_gradients |
                   update_quadrature_points |
                   update_JxW_values);
  // face finite element object for the side of the face in current cell
  FEFaceValues<dim> fvf(*fe, qf_rule,
                        update_values | update_gradients |
                        update_quadrature_points | update_normal_vectors |
                        update_JxW_values);
  
  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q = q_rule.size();
  const unsigned int n_qf = qf_rule.size();
  
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  
  // volumetric pre-assembly matrices
  std::vector<FullMatrix<double> >
  mass_at_qp (n_q, FullMatrix<double>(dofs_per_cell, dofs_per_cell));
  std::vector<std::vector<FullMatrix<double> > >
  stiffness_at_qp (n_q, std::vector<FullMatrix<double> > (n_dir, FullMatrix<double> (dofs_per_cell, dofs_per_cell)));
  
  for (typename DoFHandler<dim>::active_cell_iterator
       cell=dof_handler.begin_active();
       cell!=dof_handler.end(); ++cell)
  {
    if (cell->is_locally_owned())
    {
      fv.reinit(cell);
      cell->get_dof_indices (local_dof_indices);
      
      for (unsigned int qi=0; qi<n_q; ++qi)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            mass_at_qp[qi](i,j) = (fv.shape_value(i,qi) *
                                   fv.shape_value(j,qi));
      
      for (unsigned int qi=0; qi<n_q; ++qi)
        for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              stiffness_at_qp[qi][i_dir](i,j) = ((fv.shape_grad(i,qi) *
                                                  omega_i[i_dir])
                                                 *
                                                 (fv.shape_grad(j,qi) *
                                                  omega_i[i_dir]));
      break;
    }// local cell
  }// cell
  
  for (unsigned int k=0; k<n_total_ho_vars; ++k)
  {
    unsigned int g = get_component_group (k);
    unsigned int i_dir = get_direction (k);
    FullMatrix<double> local_mat (dofs_per_cell, dofs_per_cell);
    
    for (typename DoFHandler<dim>::active_cell_iterator
         cell=dof_handler.begin_active();
         cell!=dof_handler.end(); ++cell)
    {
      if (cell->is_locally_owned())
      {
        fv.reinit (cell);
        cell->get_dof_indices (local_dof_indices);
        local_mat = 0;
        unsigned int material_id = cell->material_id ();
        if (k==0)
          pcout << "mat id " << material_id
          << ", sigt: " << all_sigt[material_id][g]
          << ", inv_sigt: " << all_inv_sigt[material_id][g] << std::endl;
        
        for (unsigned int qi=0; qi<n_q; ++qi)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
            {
              local_mat(i,j) += (stiffness_at_qp[qi][i_dir](i,j) *
                                 all_inv_sigt[material_id][g]
                                 +
                                 mass_at_qp[qi](i,j) *
                                 all_sigt[material_id][g]) * fv.JxW(qi);
            }
        
        if (cell->at_boundary())
          for (unsigned int fn=0; fn<GeometryInfo<dim>::faces_per_cell; ++fn)
          {
            if (cell->at_boundary(fn) &&
                !is_reflective_bc[cell->face(fn)->boundary_id()])
            {
              fvf.reinit (cell,fn);
              const Tensor<1,dim> vec_n = fvf.normal_vector(0);
              double absndo = std::fabs (vec_n * omega_i[i_dir]);
              for (unsigned int qi=0; qi<n_qf; ++qi)
                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  for (unsigned int j=0; j<dofs_per_cell; ++j)
                    local_mat(i,j) += (absndo *
                                       fvf.shape_value(i,qi) *
                                       fvf.shape_value(j,qi) *
                                       fvf.JxW(qi));
            }// non-ref bd
            else if (cell->at_boundary(fn) &&
                     is_reflective_bc[cell->face(fn)->boundary_id()] &&
                     is_explicit_reflective)
            {
              fvf.reinit (cell,fn);
              unsigned int boundary_id = cell->face(fn)->boundary_id ();
              unsigned int r_dir = get_reflective_direction_index (boundary_id, i_dir);
              const Tensor<1, dim> vec_n = fvf.normal_vector (0);
              double ndo = omega_i[i_dir] * vec_n;
              for (unsigned int qi=0; qi<n_qf; ++qi)
                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  for (unsigned int j=0; j<dofs_per_cell; ++j)
                    local_mat(i,j) += (ndo *
                                       fvf.shape_value(i,qi) *
                                       (omega_i[r_dir] * fvf.shape_grad(j,qi)) *
                                       all_inv_sigt[material_id][g] *
                                       fvf.JxW(qi));
            }// explicit ref bd face
          }
        
        vec_ho_sys[k]->add (local_dof_indices,
                            local_dof_indices,
                            local_mat);
      }// local cell
    }// cell
    
    vec_ho_sys[k]->compress (VectorOperation::add);
  }// components
}

template <int dim>
void EP_SN<dim>::local_matrix_check (FullMatrix<double> &local_mat,
                                     std::string str,
                                     unsigned int ind)
{
  std::cout << str << ", " << local_mat.l1_norm() << ", ind " << ind << std::endl;
}

template <int dim>
void EP_SN<dim>::local_radio (std::string str)
{
  std::cout << str << " on proc "
  << Utilities::MPI::this_mpi_process (mpi_communicator) << std::endl;
}

template <int dim>
void EP_SN<dim>::local_radio (std::string str,
                              double &num)
{
  pcout << str << ": " << num << std::endl;
}

template <int dim>
void EP_SN<dim>::local_radio (std::string str,
                              unsigned int &num)
{
  pcout << str << ": " << num << std::endl;
}

template <int dim>
void EP_SN<dim>::assemble_ho_interface_new ()
{
  const QGauss<dim-1>  qf_rule(p_order+1);
  
  // face finite element object for the side of the face in current cell
  FEFaceValues<dim> fvf(*fe, qf_rule,
                        update_values | update_gradients |
                        update_quadrature_points | update_normal_vectors |
                        update_JxW_values);
  // face finite element object for the side of the face in neighbor cell
  FEFaceValues<dim> fvf_nei(*fe, qf_rule,
                            update_values | update_gradients |
                            update_quadrature_points | update_normal_vectors |
                            update_JxW_values);
  
  
  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_qf = qf_rule.size();
  
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  std::vector<types::global_dof_index> neigh_dof_indices (dofs_per_cell);
  
  FullMatrix<double> vp_up (dofs_per_cell, dofs_per_cell);
  FullMatrix<double> vp_un (dofs_per_cell, dofs_per_cell);
  FullMatrix<double> vn_up (dofs_per_cell, dofs_per_cell);
  FullMatrix<double> vn_un (dofs_per_cell, dofs_per_cell);
  
  for (unsigned int k=0; k<n_total_ho_vars; ++k)
  {
    unsigned int g = get_component_group (k);
    unsigned int i_dir = get_direction (k);
    
    for (typename DoFHandler<dim>::active_cell_iterator
         cell=dof_handler.begin_active();
         cell!=dof_handler.end(); ++cell)
      if (cell->is_locally_owned())
      {
        unsigned int material_id = cell->material_id ();
        double local_sigt = all_sigt[material_id][g];
        double local_inv_sigt = all_inv_sigt[material_id][g];
        double local_measure = cell->measure ();
        cell->get_dof_indices (local_dof_indices);
        
        for (unsigned int fn=0; fn<GeometryInfo<dim>::faces_per_cell; ++fn)
        {
          if (!cell->at_boundary(fn) &&
              cell->neighbor(fn)->id()<cell->id())
          {
            vp_up = 0;
            vp_un = 0;
            vn_up = 0;
            vn_un = 0;
            
            fvf.reinit (cell,fn);
            const Tensor<1,dim> vec_n = fvf.normal_vector (0);
            
            typename DoFHandler<dim>::cell_iterator neig = cell->neighbor(fn);
            fvf_nei.reinit (neig, cell->neighbor_face_no(fn));
            neig->get_dof_indices (neigh_dof_indices);
            double neigh_sigt = all_sigt[neig->material_id()][g];
            double neigh_inv_sigt = all_inv_sigt[neig->material_id()][g];
            double neigh_measure = neig->measure ();
            
            double face_measure = cell->face(fn)->measure ();
            
            double avg_mfp_inv = 0.5 * (face_measure / (local_sigt * local_measure)
                                        + face_measure / (neigh_sigt * neigh_measure));
            double sige = std::max(2.5, tensor_norms[i_dir] * c_penalty * avg_mfp_inv);
            for (unsigned int qi=0; qi<n_qf; ++qi)
              for (unsigned int i=0; i<dofs_per_cell; ++i)
                for (unsigned int j=0; j<dofs_per_cell; ++j)
                {
                  double theta = 1.0;
                  vp_up(i,j) += (sige *
                                 fvf.shape_value(i,qi) *
                                 fvf.shape_value(j,qi)
                                 + theta *
                                 (-
                                  local_inv_sigt * vec_n * omega_i[i_dir] *
                                  (omega_i[i_dir] * fvf.shape_grad(i,qi)) *
                                  fvf.shape_value(j,qi) * 0.5
                                  -
                                  local_inv_sigt * vec_n * omega_i[i_dir] *
                                  fvf.shape_value(i,qi) * 0.5 *
                                  (omega_i[i_dir] * fvf.shape_grad(j,qi)))
                                 ) * fvf.JxW(qi);
                  
                  vp_un(i,j) += (-sige *
                                 fvf.shape_value(i,qi) *
                                 fvf_nei.shape_value(j,qi)
                                 + theta *
                                 (+
                                  local_inv_sigt * vec_n * omega_i[i_dir] *
                                  (omega_i[i_dir] * fvf.shape_grad(i,qi)) *
                                  fvf_nei.shape_value(j,qi) * 0.5
                                  -
                                  neigh_inv_sigt * vec_n * omega_i[i_dir] *
                                  fvf.shape_value(i,qi) * 0.5 *
                                  (omega_i[i_dir] * fvf_nei.shape_grad(j,qi)))
                                 ) * fvf.JxW(qi);
                  
                  vn_up(i,j) += (-sige *
                                 fvf_nei.shape_value(i,qi) *
                                 fvf.shape_value(j,qi)
                                 + theta *
                                 (-
                                  neigh_inv_sigt * vec_n * omega_i[i_dir] *
                                  (omega_i[i_dir] * fvf_nei.shape_grad(i,qi)) *
                                  fvf.shape_value(j,qi) * 0.5
                                  +
                                  local_inv_sigt * vec_n * omega_i[i_dir] *
                                  fvf_nei.shape_value(i,qi) * 0.5 *
                                  (omega_i[i_dir] * fvf.shape_grad(j,qi)))
                                 ) * fvf.JxW(qi);
                  
                  vn_un(i,j) += (sige *
                                 fvf_nei.shape_value(i,qi) *
                                 fvf_nei.shape_value(j,qi)
                                 + theta *
                                 (+
                                  neigh_inv_sigt * vec_n * omega_i[i_dir] *
                                  (omega_i[i_dir] * fvf_nei.shape_grad(i,qi)) *
                                  fvf_nei.shape_value(j,qi) * 0.5
                                  +
                                  neigh_inv_sigt * vec_n * omega_i[i_dir] *
                                  fvf_nei.shape_value(i,qi) * 0.5 *
                                  (omega_i[i_dir] * fvf_nei.shape_grad(j,qi)))
                                 ) * fvf.JxW(qi);
                }
            
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
        }// face
      }// local cell
    vec_ho_sys[k]->compress(VectorOperation::add);
  }// component
}

template <int dim>
void EP_SN<dim>::initialize_ho_preconditioners ()
{
  //TimerOutput::Scope t (computing_timer, "HO preconditioner initialization");
  pcout << "tot vars " << n_total_ho_vars << std::endl;
  pre_ho_amg.resize (n_total_ho_vars);
  for (unsigned int i=0; i<n_total_ho_vars; ++i)
  {
    pre_ho_amg[i].reset ();
    pre_ho_amg[i] = (std_cxx11::shared_ptr<LA::MPI::PreconditionAMG> (new LA::MPI::PreconditionAMG));
    LA::MPI::PreconditionAMG::AdditionalData data;
    if (have_reflective_bc && is_explicit_reflective)
      data.symmetric_operator = false;
    else
      data.symmetric_operator = true;
    pre_ho_amg[i]->initialize(*(vec_ho_sys)[i], data);
  }
}

template <int dim>
void EP_SN<dim>::ho_solve ()
{
  //TimerOutput::Scope t(computing_timer, "HO solve");
  
  for (unsigned int i=0; i<n_total_ho_vars; ++i)
  {
    SolverControl solver_control (dof_handler.n_dofs(),
                                  1.0e-15);
    //vec_ho_rhs[i]->l1_norm()*1.0e-15);
    PETScWrappers::PreconditionNone precond (*(vec_ho_sys)[i]);
    
    //#ifdef USE_PETSC_LA
    pcout << "solvers " << std::endl;
    if (have_reflective_bc && is_explicit_reflective)
    {
      PETScWrappers::SolverBicgstab
      solver (solver_control, mpi_communicator);
      solver.solve (*(vec_ho_sys)[i],
                    *(vec_aflx)[i],
                    *(vec_ho_rhs)[i],
                    *(pre_ho_amg)[i]);
    }
    else
    {
      
      pcout << "solver cg" << std::endl;
      LA::SolverCG solver (solver_control, mpi_communicator);
      *(vec_aflx)[i] = 0;
      solver.solve (*(vec_ho_sys)[i],
                    *(vec_aflx)[i],
                    *(vec_ho_rhs)[i],
                    *(pre_ho_amg)[i]);
      pcout << "   Solved in " << solver_control.last_step() << std::endl;
    }
    //#else
    // LA::SolverCG solver(solver_control);
    //#endif
    //pcout << "sys norm dir " << i << ": " << vec_ho_sys[i]->l1_norm () << std::endl;
    //pcout << "rhs norm dir " << i << ": " << vec_ho_rhs[i]->l1_norm () << std::endl;
    //pcout << "aflx norm dir " << i << ": " << vec_aflx[i]->l1_norm () << std::endl;
  }
}

template <int dim>
void EP_SN<dim>::generate_moments ()
{
  // FitIt: only scalar flux is generated for now
  AssertThrow(do_nda==false, ExcMessage("Moments are generated only without NDA"));
  if (!do_nda)
    for (unsigned int g=0; g<n_group; ++g)
    {
      *(vec_ho_sflx_old)[g] = *(vec_ho_sflx)[g];
      *vec_ho_sflx[g] = 0;
      for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
        vec_ho_sflx[g]->add(wi[i_dir], *(vec_aflx)[get_component_index(i_dir, g)]);
      sflx_this_processor[g] = *vec_ho_sflx[g];
    }
}

template <int dim>
void EP_SN<dim>::generate_ho_source_new ()
{
  const QGauss<dim>  q_rule(p_order+1);
  const QGauss<dim-1>  qf_rule(p_order+1);
  
  unsigned int n_q = q_rule.size();
  unsigned int n_qf = qf_rule.size();
  
  // cell finite element object
  FEValues<dim> fv(*fe, q_rule,
                   update_values | update_gradients |
                   update_quadrature_points |
                   update_JxW_values);
  // face finite element object for the side of the face in current cell
  FEFaceValues<dim> fvf(*fe, qf_rule,
                        update_values | update_gradients |
                        update_quadrature_points | update_normal_vectors |
                        update_JxW_values);
  
  // cell rhs's
  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  Vector<double> cell_rhs (dofs_per_cell);
  
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  for (unsigned int k=0; k<n_total_ho_vars; ++k)
  {
    unsigned int g = get_component_group (k);
    unsigned int i_dir = get_direction (k);
    *(vec_ho_rhs)[k] = *(vec_ho_fixed_rhs)[g];
    
    for (typename DoFHandler<dim>::active_cell_iterator
         cell = dof_handler.begin_active();
         cell!= dof_handler.end(); ++cell)
      if (cell->is_locally_owned())
      {
        fv.reinit (cell);
        cell->get_dof_indices (local_dof_indices);
        unsigned int material_id = cell->material_id ();
        cell_rhs = 0;
        
        std::vector<std::vector<double> >
        sflx_at_qp (n_group, std::vector<double> (n_q));
        for (unsigned int gin=0; gin<n_group; ++gin)
          fv.get_function_values (sflx_this_processor[gin], sflx_at_qp[gin]);
        
        for (unsigned int qi=0; qi<n_q; ++qi)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            double test_func_jxw = (fv.shape_value(i,qi) *
                                    fv.JxW(qi));
            for (unsigned int gin=0; gin<n_group; ++gin)
              cell_rhs(i) += (test_func_jxw *
                              all_sigs_per_ster[material_id][gin][g]*
                              sflx_at_qp[gin][qi]);
          }
        
        if (have_reflective_bc && !is_explicit_reflective)
          for (unsigned int fn=0; fn<GeometryInfo<dim>::faces_per_cell; ++fn)
            if (cell->at_boundary(fn) &&
                is_reflective_bc[cell->face(fn)->boundary_id()])
            {
              local_radio("sth wrong");
              fvf.reinit (cell,fn);
              unsigned int boundary_id = cell->face(fn)->boundary_id ();
              const Tensor<1,dim> vec_n = fvf.normal_vector (0);
              unsigned int r_dir = get_reflective_direction_index (boundary_id, i_dir);
              std::vector<Tensor<1, dim> > gradients_at_qp (n_qf);
              fvf.get_function_gradients (sflx_this_processor[g], gradients_at_qp);
              
              for (unsigned int qi=0; qi<n_qf; ++qi)
                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  cell_rhs(i) += (fvf.shape_value(i, qi) *
                                  vec_n * omega_i[i_dir] *
                                  all_inv_sigt[material_id][g] *
                                  omega_i[r_dir] * gradients_at_qp[qi] *
                                  fvf.JxW(qi));
            }
        constraints.distribute_local_to_global (cell_rhs,
                                                local_dof_indices,
                                                *vec_ho_rhs[k]);
        //vec_ho_rhs[k]->add(local_dof_indices,
        //                   cell_rhs);
      }// local cell
    
    vec_ho_rhs[k]->compress (VectorOperation::add);
  }// component
}

template <int dim>
void EP_SN<dim>::NDA_PI ()
{
}

template <int dim>
void EP_SN<dim>::NDA_SI ()
{
}

template <int dim>
void EP_SN<dim>::scale_fiss_transfer_matrices ()
{
  if (do_nda)
  {
  }
  else
  {
    ho_scaled_fiss_transfer_per_ster.resize (n_material);
    for (unsigned int m=0; m<n_material; ++m)
    {
      std::vector<std::vector<double> >  tmp (n_group, std::vector<double>(n_group));
      if (is_material_fissile[m])
        for (unsigned int gin=0; gin<n_group; ++gin)
          for (unsigned int g=0; g<n_group; ++g)
            tmp[gin][g] = all_ksi_nusigf_per_ster[m][gin][g] / k_ho;
      ho_scaled_fiss_transfer_per_ster[m] = tmp;
    }
  }
}

template <int dim>
void EP_SN<dim>::generate_fixed_source_new ()
{
  const QGauss<dim>  q_rule(p_order+1);
  
  unsigned int n_q = q_rule.size();
  
  // cell finite element object
  FEValues<dim> fv(*fe, q_rule,
                   update_values | update_gradients |
                   update_quadrature_points |
                   update_JxW_values);
  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  // cell rhs's
  Vector<double> cell_rhs (dofs_per_cell);
  
  for (unsigned int g=0; g<n_group; ++g)
  {
    for (typename DoFHandler<dim>::active_cell_iterator
         cell = dof_handler.begin_active();
         cell!= dof_handler.end(); ++cell)
      if (cell->is_locally_owned())
      {
        int material_id = cell->material_id ();
        if (is_eigen_problem)
        {
          if (is_material_fissile[material_id])
          {
            fv.reinit (cell);
            cell_rhs = 0;
            cell->get_dof_indices (local_dof_indices);
            std::vector<std::vector<double> >
            local_ho_sflxes (n_group, std::vector<double> (n_q));
            
            for (unsigned int gin=0; gin<n_group; ++gin)
              fv.get_function_values (*(vec_ho_sflx)[gin], local_ho_sflxes[gin]);
            
            for (unsigned int qi=0; qi<n_q; ++qi)
              for (unsigned int i=0; i<dofs_per_cell; ++i)
              {
                double test_func_jxw = fv.shape_value (i, qi) * fv.JxW (qi);
                for (unsigned int gin=0; g<n_group; ++gin)
                  cell_rhs(i) += (test_func_jxw *
                                  ho_scaled_fiss_transfer_per_ster[material_id][gin][g] *
                                  local_ho_sflxes[gin][qi]);
              }
            
            vec_ho_fixed_rhs[g]->add(local_dof_indices,
                                     cell_rhs);
          }// fissile materials
        }// eigenvalue problem
        else
        {
          auto it = std::max_element (all_q_per_ster[material_id].begin(),
                                      all_q_per_ster[material_id].end());
          if (*it>1.0e-13)
          {
            fv.reinit (cell);
            cell->get_dof_indices (local_dof_indices);
            cell_rhs = 0.0;
            for (unsigned int qi=0; qi<n_q; ++qi)
            {
              double test_func_jxw;
              for (unsigned int i=0; i<dofs_per_cell; ++i)
              {
                test_func_jxw = fv.shape_value (i, qi) * fv.JxW (qi);
                if (all_q_per_ster[material_id][g]>1.0e-13)
                  cell_rhs(i) += test_func_jxw * all_q_per_ster[material_id][g];
              }
            }
            
            if (do_nda)
            {
              AssertThrow (!do_nda,
                           ExcMessage("NDA is for future project"));
            }
            vec_ho_fixed_rhs[g]->add(local_dof_indices,
                                     cell_rhs);
          }// nonzero source region
        }// non-eigen problem
      }// local cells
    
    vec_ho_fixed_rhs[g]->compress (VectorOperation::add);
    if (do_nda)
      vec_lo_fixed_rhs[g]->compress (VectorOperation::add);
  }// component
}

template <int dim>
void EP_SN<dim>::power_iteration ()
{
  k_ho = 1.0;
  double err_k = 1.0;
  double err_phi = 1.0;
  
  initialize_ho_preconditioners ();
  
  while (err_k>err_k_tol && err_phi>err_phi_tol)
  {
    k_ho_prev_gen = k_ho;
    
    for (unsigned int g=0; g<n_group; ++g)
      *(vec_ho_sflx_prev_gen)[g] = *(vec_ho_sflx)[g];
    
    source_iteration ();
    
    fission_source_prev_gen = fission_source;
    
    fission_source = estimate_fiss_source (vec_ho_sflx);
    
    k_ho = estimate_k (fission_source, fission_source_prev_gen, k_ho_prev_gen);
    
    double norm_factor = vec_ho_sflx[0]->l1_norm ();
    renormalize_sflx (vec_ho_sflx, norm_factor);
    
    err_phi = estimate_phi_diff (vec_ho_sflx, vec_ho_sflx_prev_gen);
    
    err_k = std::fabs (k_ho - k_ho_prev_gen) / k_ho;
  }
}

template <int dim>
void EP_SN<dim>::source_iteration ()
{
  unsigned int ct = 0;
  double err_phi = 1.0;
  double err_phi_old;
  generate_moments ();
  while (err_phi>1.0e-7)
  {
    //generate_ho_source ();
    ct += 1;
    generate_ho_source_new ();
    ho_solve ();
    generate_moments ();
    err_phi_old = err_phi;
    err_phi = estimate_phi_diff (vec_ho_sflx, vec_ho_sflx_old);
    local_radio ("iteration", ct);
    local_radio ("SI phi err", err_phi);
    double spectral_radius = err_phi / err_phi_old;
    local_radio ("spectral radius", spectral_radius);
    postprocess ();
  }
}

template <int dim>
void EP_SN<dim>::renormalize_sflx (std::vector<LA::MPI::Vector*> &target_sflxes, double &normalization_factor)
{
  AssertThrow (target_sflxes.size()==n_group,
               ExcMessage("vector of scalar fluxes must have a size of n_group"));
  for (unsigned int g=0; g<n_group; ++g)
    *(target_sflxes)[g] /= normalization_factor;
}

template <int dim>
void EP_SN<dim>::postprocess ()
{
  std::vector<double> r_collision (n_group, 0.0);
  const QGauss<dim>  q_rule(p_order+1);
  unsigned int n_q = q_rule.size ();
  FEValues<dim> fv(*fe, q_rule,
                   update_values |
                   update_quadrature_points |
                   update_JxW_values);
  
  for (typename DoFHandler<dim>::active_cell_iterator
       cell=dof_handler.begin_active();
       cell!=dof_handler.end(); ++cell)
    if (cell->is_locally_owned() &&
        cell->material_id()==0)
    {
      fv.reinit(cell);
      for (unsigned int g=0; g<n_group; ++g)
      {
        std::vector<double> local_sflx (n_q);
        fv.get_function_values (sflx_this_processor[g], local_sflx);
        for (unsigned int qi=0; qi<n_q; ++qi)
          r_collision[g] += (all_sigt[cell->material_id()][g] *
                             local_sflx[qi] *
                             fv.JxW(qi));
      }
    }
  
  std::vector<double> r_glob (n_group, 0.0);
  for (unsigned int g=0; g<n_group; ++g)
    r_glob[g] = Utilities::MPI::sum (r_collision[g], mpi_communicator);
  
  
  
  if (Utilities::MPI::this_mpi_process(mpi_communicator)==0)
  {
    std::ofstream pp;
    std::ostringstream os;
    os << discretization << "-info-" << global_refinements << ".txt";
    pp.open(os.str());
    pp << "number of cells: " << triangulation.n_global_active_cells() << std::endl;
    pp << "dof counts: " << dof_handler.n_dofs () << std::endl;
    pp << "Group | collision rates: " << std::endl;
    for (unsigned int g=0; g<n_group; ++g)
      pp << g + 1 << "    | " << std::setprecision (15) << r_glob[g] << std::endl;
    pp.close ();
  }
}

template <int dim>
double EP_SN<dim>::estimate_fiss_source (std::vector<LA::MPI::Vector*> &phis)
{
  double fiss_source = 0.0;
  
  const QGauss<dim>  q_rule(p_order+1);
  unsigned int n_q = q_rule.size ();
  FEValues<dim> fv(*fe, q_rule,
                   update_values |
                   update_quadrature_points |
                   update_JxW_values);
  
  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin (),
  endc = dof_handler.end ();
  for (; cell!=endc; ++cell)
    if (cell->is_locally_owned() &&
        is_material_fissile[cell->material_id ()])
    {
      fv.reinit (cell);
      std::vector<std::vector<double> > local_phis (n_group,
                                                    std::vector<double> (n_q));
      unsigned int material_id = cell->material_id ();
      for (unsigned int g=0; g<n_group; ++g)
        fv.get_function_values (*(phis)[g], local_phis[g]);
      
      for (unsigned int qi=0; qi<n_q; ++qi)
        for (unsigned int g=0; g<n_group; ++g)
          fiss_source += (all_nusigf[material_id][g] *
                          local_phis[g][qi] *
                          fv.JxW(qi));
    }
  // broadcasting to all processors
  return Utilities::MPI::sum (fiss_source, mpi_communicator);
}

template <int dim>
double EP_SN<dim>::estimate_k (double &fiss_source,
                               double &fiss_source_prev_gen,
                               double &k_prev_gen)
{
  // do we have to re-normalize the scalar fluxes?
  return k_prev_gen * fiss_source_prev_gen / fiss_source;
}

template <int dim>
double EP_SN<dim>::estimate_phi_diff (std::vector<LA::MPI::Vector*> &phis_newer,
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
void EP_SN<dim>::do_iterations ()
{
  initialize_ho_preconditioners ();
  generate_fixed_source_new ();
  
  if (is_eigen_problem)
  {
    if (do_nda)
      NDA_PI ();
    else
      power_iteration ();
  }
  else
  {
    if (do_nda)
      NDA_SI ();
    else
      source_iteration ();
  }
}

template <int dim>
void EP_SN<dim>::output_results () const
{
  std::string sec_name = "Graphical output";
  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  
  for (unsigned int g=0; g<n_group; ++g)
  {
    std::ostringstream os;
    os << "ho_phi_g_" << g;
    data_out.add_data_vector (sflx_this_processor[g], os.str ());
  }
  
  Vector<float> subdomain (triangulation.n_active_cells ());
  for (unsigned int i=0; i<subdomain.size(); ++i)
    subdomain(i) = triangulation.locally_owned_subdomain ();
  data_out.add_data_vector (subdomain, "subdomain");
  
  data_out.build_patches ();
  
  const std::string filename = ("sflx-" + Utilities::int_to_string
                                (triangulation.locally_owned_subdomain (), 4));
  std::ofstream output ((filename + ".vtu").c_str ());
  data_out.write_vtu (output);
  
  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
  {
    std::vector<std::string> filenames;
    for (unsigned int i=0;
         i<Utilities::MPI::n_mpi_processes(mpi_communicator);
         ++i)
      filenames.push_back ("sflx-" +
                           Utilities::int_to_string (i, 4) +
                           ".vtu");
    std::ostringstream os;
    os << "solution-" << discretization << "-" << global_refinements << ".pvtu";
    std::ofstream master_output ((os.str()).c_str ());
    data_out.write_pvtu_record (master_output, filenames);
  }
}

template <int dim>
void EP_SN<dim>::global_matrix_check (unsigned int ind)
{
  pcout << "global system norm index " << ind << ", " << vec_ho_sys[ind]->l1_norm () << std::endl;
}

template <int dim>
void EP_SN<dim>::run ()
{
  generate_globally_refined_grid ();
  setup_boundary_ids ();
  initialize_material_id ();
  setup_system ();
  report_system ();
  assemble_ho_system ();
  do_iterations ();
  
  if (Utilities::MPI::n_mpi_processes(mpi_communicator) <= 32)
  {
    pcout << "output " << std::endl;
    output_results();
  }
}

// explicit instantiation to avoid linking error
template class EP_SN<2>;
template class EP_SN<3>;
