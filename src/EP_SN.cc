/*
 namespace LA
 {
 using namespace dealii::LinearAlgebraPETSc;
 }
 */

//#include "problem_definition.h"
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
//paras(new ProblemDefinition<dim>(prm)),
err_k_tol(1.0e-6),
err_phi_tol(1.0e-6),
pcout(std::cout,
      (Utilities::MPI::this_mpi_process(mpi_communicator)
       == 0))/*,
              computing_timer (mpi_communicator,
              pcout,
              TimerOutput::summary,
              TimerOutput::wall_times)*/
{
  this->process_input ();
}

template <int dim>
EP_SN<dim>::~EP_SN()
{
  dof_handler.clear();
  delete fe;
  //delete paras;
}

template <int dim>
void EP_SN<dim>::process_input ()
{
  {
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
    all_sigs_per_ster = this->get_sigma_s ();
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
  {
    if (cell->is_locally_owned())
    {
      Point<dim> center = cell->center ();
      std::vector<unsigned int> relative_position (3);
      get_cell_relative_position (center, relative_position);
      unsigned int material_id = relative_position_to_id[relative_position];
      cell->set_material_id (material_id);
    }
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
  //TimerOutput::Scope t(computing_timer, "setup HO system");

  if (boost::iequals(discretization,"DFEM") || boost::iequals(discretization,"DG"))
    fe = new FE_DGQ<dim> (p_order);
  else
    fe = new FE_Q<dim> (p_order);

  dof_handler.distribute_dofs (*fe);

  local_dofs = dof_handler.locally_owned_dofs ();

  DynamicSparsityPattern dsp (local_dofs);

  if (discretization=="DFEM")
    DoFTools::make_flux_sparsity_pattern (dof_handler, dsp);
  else
    DoFTools::make_sparsity_pattern (dof_handler, dsp);

  // be careful with the following
  SparsityTools::distribute_sparsity_pattern (dsp,
                                              dof_handler.n_locally_owned_dofs_per_processor (),
                                              mpi_communicator,
                                              local_dofs);

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

  c_penalty = p_order * (p_order + 1.0);
}

/*
 template <int dim>
 void EP_SN<dim>::assemble_diffusion()
 {
 // TimerOutput::Scope t(computing_timer, "assembly LO");

 const QGauss<dim> q_rule (2);

 FEValues<dim> fv(*lo_fe,q_rule,
 update_values | update_gradients |
 update_quadrature_points | update_JxW_values);

 const unsigned int dofs_per_cell = lo_fe->dofs_per_cell;
 const unsigned int n_q = q_rule.size();

 FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 FullMatrix<double> cell_streaming_matrix(dofs_per_cell, dofs_per_cell);
 FullMatrix<double> cell_collision_matrix(dofs_per_cell, dofs_per_cell);

 Vector<double> cell_rhs(dofs_per_cell);

 std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

 std::vector<FullMatrix<double> > vec_collision_local (n_q, FullMatrix<double> (dofs_per_cell, dofs_per_cell));
 std::vector<FullMatrix<double> > vec_streaming_local (n_q, FullMatrix<double> (dofs_per_cell, dofs_per_cell));
 int pre_assemble_cell = 0;

 typename DoFHandler<dim>::active_cell_iterator
 cell = dof_handler.begin_active(),
 endc = dof_handler.end();
 for (; cell!=endc; ++cell,++i_cell)
 {
 if (cell->is_locally_owned())
 {
 fv.reinit(cell);
 cell_matrix = 0.0;
 cell_collision_matrix = 0.0;
 cell_streaming_matrix = 0.0;

 if (pre_assemble_cell==0)
 {
 for (unsigned int qi=0; qi<n_q; ++qi)
 for (unsigned int i=0; i<dofs_per_cell; ++i)
 for (unsigned int j=0; j<dofs_per_cell; ++j)
 vec_collision_local[qi](i,j) += (fv.shape_value(i,qi) *
 fv.shape_value(j,qi));


 for (unsigned int qi=0; qi<n_q; ++qi)
 for (unsigned int i=0; i<dofs_per_cell; ++i)
 for (unsigned int j=0; j<dofs_per_cell; ++j)
 vec_streaming_local[qi][g](i,j) += (fv.shape_grad(i,qi) *
 fv.shape_grad(j,qi));

 pre_assemble_cell = 1;
 }// assemble once on every processor for diffusion streaming and collision

 for (unsigned int g=0; g<n_group; ++g)
 {
 double cell_siga = sigt[g][cell->material_id()] - sig0[g][cell->material_id()];
 double cell_dif_coeff = cell_dif_coeff[g][cell->material_id()];

 for (unsigned int qi=0; qi<n_q; ++qi)
 cell_matrix.add(cell_dif_coeff,vec_streaming_local[qi],cell_siga,vec_collision_local[qi]);


 }// loop over all groups

 }// locally owned cells
 }// cell
 }
 */

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
  assemble_ho_volume_boundary ();

  if (discretization=="DFEM")
  {
    local_radio ("Assemble cell interface bilinear forms for DFEM");
    assemble_ho_interface ();
  }
}

template <int dim>
void EP_SN<dim>::assemble_ho_volume_boundary ()
{
  //TimerOutput::Scope t(computing_timer, "assembly HO");

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

  bool pre_assemble_cell_finished = false;

  for (typename DoFHandler<dim>::active_cell_iterator cell=dof_handler.begin_active();
       cell!=dof_handler.end(); ++cell)
  {
    if (cell->is_locally_owned())
    {
      fv.reinit(cell);
      cell->get_dof_indices (local_dof_indices);

      std::vector<double> local_sigt = all_sigt[cell->material_id()];
      std::vector<double> local_inv_sigt = all_inv_sigt[cell->material_id()];

      std::vector<FullMatrix<double> >
      local_matrices (n_total_ho_vars, FullMatrix<double> (dofs_per_cell, dofs_per_cell));
      // FixIt: a more proper definition for h_cell

      if (!pre_assemble_cell_finished)
      {
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

        pre_assemble_cell_finished = true;
      }

      // Use pre-assembled matrix components in reference cell to do assembly in real matrix
      FullMatrix<double>
      unscaled_mass (dofs_per_cell, dofs_per_cell);
      std::vector<FullMatrix<double> >
      unscaled_stiffness(n_dir, FullMatrix<double>(dofs_per_cell, dofs_per_cell));

      for (unsigned int qi=0; qi<n_q; ++qi)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            unscaled_mass(i,j) += mass_at_qp[qi](i,j) * fv.JxW(qi);

      for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
        for (unsigned int qi=0; qi<n_q; ++qi)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              unscaled_stiffness[i_dir](i,j) += (stiffness_at_qp[qi][i_dir](i,j) *
                                                 fv.JxW(qi));

      for (unsigned int g=0; g<n_group; ++g)
      {
        //FullMatrix<double> scaled_mass (dofs_per_cell, dofs_per_cell);
        FullMatrix<double> scaled_mass = unscaled_mass;
        scaled_mass *= local_sigt[g];
        for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
        {
          unsigned int ind = get_component_index (i_dir, g);
          local_matrices[ind] = scaled_mass;
          local_matrices[ind].add (local_inv_sigt[g], unscaled_stiffness[i_dir]);
        }
      }

      for (unsigned int fn=0; fn<GeometryInfo<dim>::faces_per_cell; ++fn)
      {
        if (cell->face(fn)->at_boundary())
        {
          fvf.reinit (cell, fn);
          const Tensor<1, dim> vec_n = fvf.normal_vector (0);
          unsigned int boundary_id = cell->face(fn)->boundary_id ();
          if (!is_reflective_bc[boundary_id])
          {
            std::ostringstream os;
            os << "assemble bd " << boundary_id;
            local_radio (os.str());
            for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
            {
              double absndo = std::fabs (vec_n * omega_i[i_dir]);
              // Note: we assume the face is not curvilinear
              for (unsigned int g=0; g<n_group; ++g)
              {
                int ind = get_component_index (i_dir, g);
                for (unsigned int qi=0; qi<n_qf; ++qi)
                  for (unsigned int i=0; i<dofs_per_cell; ++i)
                    for (unsigned int j=0; j<dofs_per_cell; ++j)
                      local_matrices[ind](i,j) += (absndo *
                                                   fvf.shape_value(i,qi) *
                                                   fvf.shape_value(j,qi) *
                                                   fvf.JxW(qi));


              }// g
            }// i_dir
          }// non-reflective boundary
          else
          {
            std::cout << "sth is wrong" << std::endl;
            if (is_explicit_reflective)// assemble nothing if false
              for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
              {
                unsigned int r_dir = get_reflective_direction_index (boundary_id, i_dir);
                double absndo = omega_i[i_dir] * vec_n;
                for (unsigned int g=0; g<n_group; ++g)
                {
                  unsigned int ind = get_component_index (i_dir, g);
                  for (unsigned int qi=0; qi<n_qf; ++qi)
                    for (unsigned int i=0; i<dofs_per_cell; ++i)
                      for (unsigned int j=0; j<dofs_per_cell; ++j)
                        local_matrices[ind](i,j) += (absndo *
                                                     fvf.shape_value(i,qi) *
                                                     (omega_i[r_dir] * fvf.shape_grad(j,qi)) /
                                                     local_sigt[g] *
                                                     fvf.JxW(qi));
                }// g
              }// i_dir
          }// is_reflective_bc
        }// boundary face
      }// face

      for (unsigned int k=0; k<n_total_ho_vars; ++k)
      {
        std::ostringstream os;
        os << "loc_to_glob mapping cell " << cell->id() << " component " << k;
        local_radio (os.str());
        vec_ho_sys[k]->add (local_dof_indices,
                            local_dof_indices,
                            local_matrices[k]);
        //pcout << "local_matrices " << k << " norm " << local_matrices[k].l1_norm() << std::endl;
        //pcout << "sys norms: " << vec_ho_sys[k]->l1_norm () << " on proc "
        //<< Utilities::MPI::this_mpi_process(mpi_communicator) << std::endl;
      }

      // FixIt: constraints are needed if refinement work is desired
    }// cell locally owned
  }// cell

  for (unsigned int k=0; k<1; ++k)
  {
    std::ostringstream os;
    os << "compress sys_mats begins component " << k;
    local_radio (os.str());
    vec_ho_sys[k]->compress (VectorOperation::add);
    std::ostringstream os1;
    os1 << "compress sys_mats ends component " << k;
    local_radio (os1.str());

  }
  for (unsigned int k=1; k<2; ++k)
  {
    std::ostringstream os;
    os << "compress sys_mats begins component " << k;
    local_radio (os.str());
    vec_ho_sys[k]->compress (VectorOperation::add);
    std::ostringstream os1;
    os1 << "compress sys_mats ends component " << k;
    local_radio (os1.str());

  }
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
void EP_SN<dim>::assemble_ho_interface ()
{
  //TimerOutput::Scope t(computing_timer, "assembly HO");

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

  // face terms: v^\pm * u^\pm
  std::vector<FullMatrix<double> > all_real_vp_up(n_group * n_dir, FullMatrix<double>(dofs_per_cell, dofs_per_cell));
  std::vector<FullMatrix<double> > all_real_vp_un(n_group * n_dir, FullMatrix<double>(dofs_per_cell, dofs_per_cell));
  std::vector<FullMatrix<double> > all_real_vn_up(n_group * n_dir, FullMatrix<double>(dofs_per_cell, dofs_per_cell));
  std::vector<FullMatrix<double> > all_real_vn_un(n_group * n_dir, FullMatrix<double>(dofs_per_cell, dofs_per_cell));

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  std::vector<types::global_dof_index> neigh_dof_indices (dofs_per_cell);

  // face pre-assembly matrices: value penalty
  std::vector<FullMatrix<double> > vec_vp_up (n_qf,
                                              FullMatrix<double> (dofs_per_cell, dofs_per_cell));
  std::vector<FullMatrix<double> > vec_vp_un (n_qf,
                                              FullMatrix<double> (dofs_per_cell, dofs_per_cell));
  std::vector<FullMatrix<double> > vec_vn_up (n_qf,
                                              FullMatrix<double> (dofs_per_cell, dofs_per_cell));
  std::vector<FullMatrix<double> > vec_vn_un (n_qf,
                                              FullMatrix<double> (dofs_per_cell, dofs_per_cell));

  // face pre-assembly matrices: gradient penalty 1
  std::vector<std::vector<FullMatrix<double> > > vec_dvp_up (n_qf,
                                                             std::vector<FullMatrix<double> > (n_dir, FullMatrix<double> (dofs_per_cell, dofs_per_cell)));
  std::vector<std::vector<FullMatrix<double> > > vec_dvp_un (n_qf,
                                                             std::vector<FullMatrix<double> > (n_dir, FullMatrix<double> (dofs_per_cell, dofs_per_cell)));
  std::vector<std::vector<FullMatrix<double> > > vec_dvn_up (n_qf,
                                                             std::vector<FullMatrix<double> > (n_dir, FullMatrix<double> (dofs_per_cell, dofs_per_cell)));
  std::vector<std::vector<FullMatrix<double> > > vec_dvn_un (n_qf,
                                                             std::vector<FullMatrix<double> > (n_dir, FullMatrix<double> (dofs_per_cell, dofs_per_cell)));

  // face pre-assembly matrices: gradient penalty 2
  std::vector<std::vector<FullMatrix<double> > > vec_vp_dup (n_qf,
                                                             std::vector<FullMatrix<double> > (n_dir, FullMatrix<double> (dofs_per_cell, dofs_per_cell)));
  std::vector<std::vector<FullMatrix<double> > > vec_vp_dun (n_qf,
                                                             std::vector<FullMatrix<double> > (n_dir, FullMatrix<double> (dofs_per_cell, dofs_per_cell)));
  std::vector<std::vector<FullMatrix<double> > > vec_vn_dup (n_qf,
                                                             std::vector<FullMatrix<double> > (n_dir, FullMatrix<double> (dofs_per_cell, dofs_per_cell)));
  std::vector<std::vector<FullMatrix<double> > > vec_vn_dun (n_qf,
                                                             std::vector<FullMatrix<double> > (n_dir, FullMatrix<double> (dofs_per_cell, dofs_per_cell)));

  bool pre_assemble_face_finished = false;

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (typename DoFHandler<dim>::active_cell_iterator cell=dof_handler.begin_active();
       cell!=dof_handler.end(); ++cell)
  {
    if (cell->is_locally_owned())
    {
      cell->get_dof_indices (local_dof_indices);
      std::vector<double> local_sigt = all_sigt[cell->material_id()];
      std::vector<double> local_inv_sigt = all_inv_sigt[cell->material_id()];
      std::vector<double> local_mfps (n_group);
      unsigned int material_id = cell->material_id ();
      double h = cell->diameter () / std::sqrt (2.0);
      get_cell_mfps (material_id, h, local_mfps);
      // FixIt: a more proper definition for h_cell

      for (unsigned int fn=0; fn<GeometryInfo<dim>::faces_per_cell; ++fn)
      {
        if (!cell->face(fn)->at_boundary() &&
            cell->neighbor(fn)->id()<cell->id())
        {
          typename DoFHandler<dim>::cell_iterator neigh = cell->neighbor(fn);
          // initialize the elements of sides of the face in current cell and neighbor
          fvf.reinit(cell, fn);
          fvf_nei.reinit(neigh, cell->neighbor_face_no(fn));
          Tensor<1,dim> n_vec = fvf.normal_vector (0);
          double sige;
          std::vector<double> neigh_sigt = all_sigt[neigh->material_id()];
          std::vector<double> neigh_mfps;
          unsigned int neigh_id = neigh->material_id ();
          double neigh_h = neigh->diameter () / std::sqrt (2.0);
          get_cell_mfps (neigh_id, neigh_h, neigh_mfps);

          if (!pre_assemble_face_finished)
          {
            // assemble once for vp/n, up/n in a reference cell
            for (unsigned int qi=0; qi<n_qf; ++qi)
            {
              double jxw = fvf.JxW(qi);
              for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
                for (unsigned int g=0; g<n_group; ++g)
                {
                  int ind = get_component_index (i_dir, g);
                  for (unsigned int i=0; i<dofs_per_cell; ++i)
                    for (unsigned int j=0; j<dofs_per_cell; ++j)
                    {
                      // ([v],\sigma_e [u])
                      if (g==0 && i_dir==0)
                      {
                        // vp_up
                        vec_vp_up[qi](i,j) += (fvf.shape_value(i,qi) *
                                               fvf.shape_value(j,qi));
                        // vp_un
                        vec_vp_un[qi](i,j) += (fvf.shape_value(i,qi) *
                                               fvf_nei.shape_value(j,qi));
                        // vn_up
                        vec_vn_up[qi](i,j) += (fvf_nei.shape_value(i,qi) *
                                               fvf.shape_value(j,qi));
                        // vn_un
                        vec_vn_un[qi](i,j) += (fvf_nei.shape_value(i,qi) *
                                               fvf_nei.shape_value(j,qi));
                      }

                      // ([v],{n*Omega*1/\sigma_t*Omega*du})
                      if (g==0)
                      {
                        // vp_dup
                        vec_vp_dup[qi][i_dir](i,j) += (fvf.shape_value(i,qi) *
                                                       0.5 * (omega_i[i_dir] * n_vec) *
                                                       (omega_i[i_dir] * fvf.shape_grad(j,qi)));
                        // vp_dun
                        vec_vp_dun[qi][i_dir](i,j) += (fvf.shape_value(i,qi) *
                                                       0.5 * (omega_i[i_dir] * n_vec) *
                                                       (omega_i[i_dir] * fvf_nei.shape_grad(j,qi)));
                        // vn_dup
                        vec_vn_dup[qi][i_dir](i,j) += (fvf_nei.shape_value(i,qi) *
                                                       0.5 * (omega_i[i_dir] * n_vec) *
                                                       (omega_i[i_dir] * fvf.shape_grad(j,qi)));
                        // vn_dun
                        vec_vn_dun[qi][i_dir](i,j) += (fvf_nei.shape_value(i,qi) *
                                                       0.5 * (omega_i[i_dir] * n_vec) *
                                                       (omega_i[i_dir] * fvf_nei.shape_grad(j,qi)));

                        // ({n*Omega*1/\sigma_t*Omega*grad_v},[u])
                        // dvp_up
                        vec_dvp_up[qi][i_dir](i,j) += (omega_i[i_dir] * fvf.shape_grad(i,qi) *
                                                       0.5 * (omega_i[i_dir] * n_vec) *
                                                       fvf.shape_value(j,qi));
                        // dvp_un
                        vec_dvp_un[qi][i_dir](i,j) += (omega_i[i_dir] * fvf.shape_grad(i,qi) *
                                                       0.5 * (omega_i[i_dir] * n_vec) *
                                                       fvf_nei.shape_value(j,qi));
                        // dvn_up
                        vec_dvn_up[qi][i_dir](i,j) += (omega_i[i_dir] * fvf_nei.shape_grad(i,qi) *
                                                       0.5 * (omega_i[i_dir] * n_vec) *
                                                       fvf.shape_value(j,qi));
                        // dvn_un
                        vec_dvn_un[qi][i_dir](i,j) += (omega_i[i_dir] * fvf_nei.shape_grad(i,qi) *
                                                       0.5 * (omega_i[i_dir] * n_vec) *
                                                       fvf_nei.shape_value(j,qi));
                      }
                    }// j
                }// g
            }// qi

            pre_assemble_face_finished = true;
          }// pre_assemble

          // Initialize all face matrices in real cells
          for (unsigned int k=0; k<n_total_ho_vars; ++k)
          {
            all_real_vp_up[k] = 0;
            all_real_vp_un[k] = 0;
            all_real_vn_up[k] = 0;
            all_real_vn_un[k] = 0;
          }
          // FixIt: try different penalty number sige

          for (unsigned int qi=0; qi<n_qf; ++qi)
          {
            double jxw = fvf.JxW(qi);
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              for (unsigned int j=0; j<dofs_per_cell; ++j)
              {
                double vp_up_jxw;
                double vp_un_jxw;
                double vn_up_jxw;
                double vn_un_jxw;
                double vp_dup_jxw;
                double vp_dun_jxw;
                double vn_dup_jxw;
                double vn_dun_jxw;
                double dvp_up_jxw;
                double dvp_un_jxw;
                double dvn_up_jxw;
                double dvn_un_jxw;

                for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
                {
                  for (unsigned int g=0; g<n_group; ++g)
                  {
                    int ind = get_component_index (i_dir, g);
                    // get cross sections ones per cell/neighbor
                    if (qi==0 && i==0 && j==0)
                      sige = std::max(0.25, (tensor_norms[i_dir] / local_mfps[g] + tensor_norms[i_dir] / neigh_mfps[g]));

                    // The following is calculating terms from reference cell to quadrature points
                    // value jump for only one group one direction
                    if (g==0 && i_dir==0)
                    {
                      vp_up_jxw = vec_vp_up[qi](i,j) * jxw;
                      vp_un_jxw = vec_vp_un[qi](i,j) * jxw;
                      vn_up_jxw = vec_vn_up[qi](i,j) * jxw;
                      vn_un_jxw = vec_vn_un[qi](i,j) * jxw;
                    }// do jxw calculation for only one group
                    if (g==0)
                    {
                      dvp_up_jxw = vec_dvp_up[qi][i_dir](i,j) * jxw;
                      dvp_un_jxw = vec_dvp_un[qi][i_dir](i,j) * jxw;
                      dvn_up_jxw = vec_dvn_up[qi][i_dir](i,j) * jxw;
                      dvn_un_jxw = vec_dvn_un[qi][i_dir](i,j) * jxw;

                      vp_dup_jxw = vec_vp_dup[qi][i_dir](i,j) * jxw;
                      vp_dun_jxw = vec_vp_dun[qi][i_dir](i,j) * jxw;
                      vn_dup_jxw = vec_vn_dup[qi][i_dir](i,j) * jxw;
                      vn_dun_jxw = vec_vn_dun[qi][i_dir](i,j) * jxw;
                    }
                    all_real_vp_up[ind](i,j) += (sige * vp_up_jxw
                                                 -
                                                 dvp_up_jxw / local_sigt[g]
                                                 -
                                                 vp_dup_jxw / local_sigt[g]);

                    all_real_vp_un[ind](i,j) += (-sige * vp_un_jxw
                                                 +
                                                 dvp_un_jxw / local_sigt[g]
                                                 -
                                                 vp_dun_jxw / neigh_sigt[g]);

                    all_real_vn_up[ind](i,j) += (-sige * vp_up_jxw
                                                 -
                                                 dvn_up_jxw / neigh_sigt[g]
                                                 +
                                                 vn_dup_jxw / local_sigt[g]);

                    all_real_vn_un[ind](i,j) += (sige * vp_up_jxw
                                                 +
                                                 dvn_un_jxw / neigh_sigt[g]
                                                 +
                                                 vn_dun_jxw / neigh_sigt[g]);
                  }// g
                }// i_dir
              }// j
          }// qi

          neigh->get_dof_indices(neigh_dof_indices);

          for (unsigned int k=0; k<n_total_ho_vars; ++k)
          {
            vec_ho_sys[k]->add (local_dof_indices,
                                local_dof_indices,
                                all_real_vp_up[k]);

            vec_ho_sys[k]->add (local_dof_indices,
                                neigh_dof_indices,
                                all_real_vp_un[k]);

            vec_ho_sys[k]->add (neigh_dof_indices,
                                local_dof_indices,
                                all_real_vn_up[k]);

            vec_ho_sys[k]->add (neigh_dof_indices,
                                neigh_dof_indices,
                                all_real_vn_un[k]);
          }
        }// non-boundary face
      }// face

      // FixIt: constraints are needed if refinement work is desired
      /*
       constraints.distribute_local_to_global (cell_matrix,
       local_dof_indices,
       system_matrix);
       */
      /*constraints.distribute_local_to_global (cell_matrix,
       cell_rhs,
       local_dof_indices,
       *(vec_sys)[0],
       system_rhs);*/
    }// cell locally owned
  }// cell

  for (unsigned int k=0; k<n_total_ho_vars; ++k)
    vec_ho_sys[k]->compress (VectorOperation::add);
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
      /*LA::SolverBicgstab solver (solver_control, mpi_communicator);
       solver.solve (*(vec_ho_sys)[i],
       *(vec_aflx)[i],
       *(vec_ho_rhs)[i],
       *(pre_ho_amg)[i]);
       */
    }
    else
    {

      pcout << "solver cg" << std::endl;
      LA::SolverCG solver (solver_control, mpi_communicator);
      *(vec_aflx)[i] = 0;
      solver.solve (*(vec_ho_sys)[i],
                    *(vec_aflx)[i],
                    *(vec_ho_rhs)[i],
                                        precond);
                    //*(pre_ho_amg)[i]);
      pcout << "   Solved in " << solver_control.last_step() << std::endl;
    }
    //#else
    // LA::SolverCG solver(solver_control);
    //#endif
    pcout << "sys norm dir " << i << ": " << vec_ho_sys[i]->l1_norm () << std::endl;
    pcout << "rhs norm dir " << i << ": " << vec_ho_rhs[i]->l1_norm () << std::endl;
    pcout << "aflx norm dir " << i << ": " << vec_aflx[i]->l1_norm () << std::endl;
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
      *(vec_ho_sflx)[g] = 0;
      *vec_ho_sflx[g]=*vec_aflx[1];
      //for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
      //  vec_ho_sflx[g]->add(wi[i_dir], *(vec_aflx)[get_component_index(i_dir, g)]);
    }
}

template <int dim>
void EP_SN<dim>::generate_ho_source ()
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
  std::vector<Vector<double> > vec_cell_rhs_reflective_bc(n_total_ho_vars, Vector<double> (dofs_per_cell));

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  for (unsigned int i=0; i<n_total_ho_vars; ++i)
    *(vec_ho_rhs)[i] = *(vec_ho_fixed_rhs)[get_component_group(i)];

  for (typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();
       cell!= dof_handler.end(); ++cell)
  {
    if (cell->is_locally_owned())
    {
      fv.reinit(cell);
      cell->get_dof_indices (local_dof_indices);
      unsigned int material_id = cell->material_id ();
      std::vector<Vector<double> > vec_cell_rhs(n_total_ho_vars,
                                                Vector<double> (dofs_per_cell));

      std::vector<std::vector<double> > all_cell_sflx(n_group, std::vector<double> (n_q));
      std::vector<double> local_sigt = all_sigt[material_id];
      for (unsigned int gin=0; gin<n_group; ++gin)
        fv.get_function_values(*(vec_ho_sflx)[gin], all_cell_sflx[gin]);

      for (unsigned int qi=0; qi<n_q; ++qi)
      {
        // do something
        double jxw = fv.JxW(qi);
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          double test_func_jxw = fv.shape_value(i,qi) * jxw;
          for (unsigned int g=0; g<n_group; ++g)
            for (unsigned int gin=0; gin<n_group; ++gin)
            {
              double addition;
              for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
              {
                unsigned int ind = get_component_index (i_dir, g);
                if (all_sigs_per_ster[material_id][gin][g]>1.0e-13 &&
                    i_dir==0)
                  addition =0;// test_func_jxw * all_sigs_per_ster[material_id][gin][g];

                vec_cell_rhs[ind](i) += addition;
              }
            }
        }
      }// qi

      if (cell->at_boundary())
      {
        // Boundary parts
        for (unsigned int fn=0; fn<GeometryInfo<dim>::faces_per_cell; ++fn)
        {
          if (cell->at_boundary (fn) &&
              (is_reflective_bc[cell->face(fn)->boundary_id()] &&
               (!is_explicit_reflective)))
          {
            unsigned int boundary_id = cell->face(fn)->boundary_id ();
            fvf.reinit(cell, fn);
            for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
            {
              unsigned int r_dir = get_reflective_direction_index (boundary_id, i_dir);
              for (unsigned int g=0; g<n_group; ++g)
              {
                int ind = get_component_index (i_dir, g);
                const Tensor<1, dim> vec_n = fvf.normal_vector (0);
                std::vector<Tensor<1, dim> > cell_daflx(n_qf);
                fvf.get_function_gradients(*(vec_aflx)[get_component_index (r_dir, g)], cell_daflx);
                for (unsigned int qi=0; qi<n_qf; ++qi)
                  for (unsigned int i=0; i<dofs_per_cell; ++i)
                    vec_cell_rhs[ind](i) += (fvf.shape_value(i, qi) *
                                             vec_n * omega_i[i_dir] / local_sigt[g] *
                                             omega_i[r_dir] * cell_daflx[qi] *
                                             fvf.JxW(qi));
              }// g
            }// i_dir
          }// reflective boundary face
        }// corresponding boundary faces
      }// cell at boundary
      for (unsigned int k=0; k<n_total_ho_vars; ++k)
        vec_ho_rhs[k]->add(local_dof_indices,
                           vec_cell_rhs[k]);

    }// local cells
  }// local

  for (unsigned int i=0; i<n_total_ho_vars; ++i)
  {
    std::ostringstream os;
    os << "sys_rhses component " << i << " compress ";
    local_radio (os.str());
    vec_ho_rhs[i]->compress (VectorOperation::add);
  }
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
void EP_SN<dim>::generate_fixed_source ()
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
  std::vector<Vector<double> > vec_cell_rhs(n_total_ho_vars,
                                            Vector<double> (dofs_per_cell));

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
          cell->get_dof_indices (local_dof_indices);
          std::vector<std::vector<double> > local_ho_sflxes (n_group, std::vector<double> (dofs_per_cell));

          for (unsigned int g=0; g<n_group; ++g)
            fv.get_function_values (*(vec_ho_sflx)[g], local_ho_sflxes[g]);

          for (unsigned int qi=0; qi<n_q; ++qi)
          {
            double test_func_jxw;
            for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              test_func_jxw = fv.shape_value (i, qi) * fv.JxW (qi);
              for (unsigned int g=0; g<n_group; ++g)
                for (unsigned int gin=0; g<n_group; ++gin)
                  vec_cell_rhs[g](i) += (test_func_jxw *
                                         ho_scaled_fiss_transfer_per_ster[material_id][gin][g]);
            }
          }

          for (unsigned int g=0; g<n_group; ++g)
            vec_ho_fixed_rhs[g]->add(local_dof_indices,
                                     vec_cell_rhs[g]);
        }
      }
      else
      {
        auto it = std::max_element (all_q_per_ster[material_id].begin(),
                                    all_q_per_ster[material_id].end());
        if (*it>1.0e-13)
        {
          fv.reinit (cell);
          cell->get_dof_indices (local_dof_indices);
          for (unsigned int qi=0; qi<n_q; ++qi)
          {
            double test_func_jxw;
            for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              test_func_jxw = fv.shape_value (i, qi) * fv.JxW (qi);
              for (unsigned int g=0; g<n_group; ++g)
                if (all_q_per_ster[material_id][g]>1.0e-13)
                  vec_cell_rhs[g](i) += test_func_jxw * all_q_per_ster[material_id][g];
            }
          }

          if (do_nda)
          {
          }
          for (unsigned int g=0; g<n_group; ++g)
            vec_ho_fixed_rhs[g]->add (local_dof_indices,
                                      vec_cell_rhs[g]);
          //for (unsigned int i=0; i<dofs_per_cell; ++i)
          //vec_ho_fixed_rhs[g]->add (local_dof_indices[i],
          //                        vec_cell_rhs[g](i));
        }
      }
    }// local cells

  for (unsigned int g=0; g<n_group; ++g)
  {
    vec_ho_fixed_rhs[g]->compress (VectorOperation::add);
    pcout << "src l1_norm " << vec_ho_fixed_rhs[g]->l1_norm () << std::endl;
    if (do_nda)
      vec_lo_fixed_rhs[g]->compress (VectorOperation::add);
  }
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
  double err_phi = 1.0;
  while (err_phi>1.0e-7)
  {
    generate_ho_source ();
    ho_solve ();
    generate_moments ();
    err_phi = estimate_phi_diff (vec_ho_sflx, vec_ho_sflx_old);
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
  generate_fixed_source ();

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
  pcout << "zero group norm " << vec_ho_sflx[0]->l1_norm () << std::endl;

  std::string sec_name = "Graphical output";
  //TimerOutput::Scope t (computing_timer, sec_name);
  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);

  for (unsigned int g=0; g<n_group; ++g)
  {
    std::ostringstream os;
    os << "ho_phi_g_" << g;
    data_out.add_data_vector (*(vec_ho_sflx)[g], os.str ());
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
    std::string base = "solution.pvtu";
    std::ofstream master_output (base.c_str ());
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
    //TimerOutput::Scope t(computing_timer, "output");
    output_results();
  }
  /*
   computing_timer.print_summary ();
   computing_timer.reset ();
   */
}

// explicit instantiation to avoid linking error
template class EP_SN<2>;
template class EP_SN<3>;
