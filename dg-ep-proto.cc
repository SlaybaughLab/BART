/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2009 - 2016 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 
 *
 * Author: Wolfgang Bangerth, Texas A&M University, 2009, 2010
 *         Timo Heister, University of Goettingen, 2009, 2010
 */

/* ---------------------------------------------------------------------
 *
 * Author: Weixiong Zheng
 *
 */
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/generic_linear_algebra.h>

namespace LA
{
  /*
   #if defined(DEAL_II_WITH_PETSC) && !(defined(DEAL_II_WITH_TRILINOS) && defined(FORCE_USE_OF_TRILINOS))
   using namespace dealii::LinearAlgebraPETSc;
   #  define USE_PETSC_LA
   #elif defined(DEAL_II_WITH_TRILINOS)
   using namespace dealii::LinearAlgebraTrilinos;
   #else
   #  error DEAL_II_WITH_PETSC or DEAL_II_WITH_TRILINOS required
   #endif
   */
  using namespace dealii::LinearAlgebraTrilinos;
  //using namespace dealii::LinearAlgebraPETSc;
}

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/cell_id.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/parameter_handler.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

#include <boost/algorithm/string.hpp>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <map>

using namespace dealii;

template <int dim>
class EP_SN
{
public:
  EP_SN (const ParameterHandler &prm, const int dimen);
  ~EP_SN();
  
  void run();
  
  static void declare_parameters (ParameterHandler &prm);
  static void declare_material_parameters (ParameterHandler &prm, ParameterHandler &prm2);
  
private:
  void initialize_index (); 
  void setup_system ();
  void generate_globally_refined_grid ();
  void report_system ();
  // void setup_lo_system();
  void setup_boundary_types ();
  void process_input_xsec ();
  void get_cell_sigts (unsigned int &material_id,
                       std::vector<double> &local_sigts);
  void get_cell_mfps (typename Triangulation<dim>::cell_iterator &cell,
                      std::vector<double> &local_sigts, std::vector<double> &local_mfps);
  void assemble_ho_system ();
  void do_iterations ();
  void assemble_ho_rhs ();
  void angular_quad ();
  void initialize_ref_bc_index ();

  unsigned int get_component_index (unsigned int &incident_angle_index, unsigned int &g);
  unsigned int get_reflective_direction_index (unsigned int &boundary_id, 
                                               unsigned int &incident_angle_index);
  
  void assemble_lo_system ();
  void prepare_correction_aflx ();
  
  void ho_solve (unsigned int &i_dir, unsigned int &g);
  void lo_solve (unsigned int g);
  void refine_grid ();
  void output_results () const;
  void power_iteration ();
  
  MPI_Comm mpi_communicator;
  
  parallel::distributed::Triangulation<dim> triangulation;
  
  /*
  DoFHandler<dim> lo_dof_handler;
  FE_DGQ<dim> *lo_fe;
  */
   
  DoFHandler<dim> dof_handler;
  FE_DGQ<dim> *fe;
  
  // std::vector<FEValuesExtractors::Scalar> comp;
  
  // FixIt: involve relevant_dofs for future if refinement is necessary
  // IndexSet lo_relevant_dofs;
  // IndexSet ho_relevant_dofs;
  IndexSet local_dofs;
  
  // ConstraintMatrix constraints;
  
  // Commented are for system-way of assembly
  /*
  LA::MPI::SparseMatrix ho_sys;
  LA::MPI::Vector ho_aflx;
  LA::MPI::Vector ho_aflx_old;
  LA::MPI::Vector delta_ho_aflx;
  LA::MPI::Vector ho_rhs;
  */
 

  // HO system
  std::vector<LA::MPI::SparseMatrix*> vec_ho_sys;
  std::vector<LA::MPI::Vector*> vec_aflx;
  std::vector<LA::MPI::Vector*> vec_ho_rhs;
  std::vector<LA::MPI::Vector*> vec_ho_sflx;
  std::vector<LA::MPI::Vector*> vec_ho_sflx_old;
  
  // LO system
  std::vector<LA::MPI::SparseMatrix*> vec_lo_sys;
  std::vector<LA::MPI::Vector*> vec_lo_rhs;
  std::vector<LA::MPI::Vector*> vec_lo_solu;
  
  std::map<std::pair<unsigned int, unsigned int>, unsigned int> component_index;
  std::map<std::pair<unsigned int, unsigned int>, unsigned int> reflective_direction_index;
  int ncells;
  unsigned int n_azi;
  unsigned int n_total_ho_vars;
  int p_order;
  double c_penalty;
  
  std::vector<Tensor<1, dim> > omega_i;
  std::vector<double> wi;
  std::vector<double> tensor_norms;
  std::vector<std::map<int, int> > ref_bc_maps;
  
  bool is_eigen_problem;
  bool do_nda;
  std::vector<bool> is_reflective_bc;
  
  ConditionalOStream                        pcout;
  TimerOutput                               computing_timer;
  
  std::vector<std_cxx11::shared_ptr<LA::MPI::PreconditionAMG> > pre_AMG;
};

template <int dim>
EP_SN<dim>::EP_SN (const ParameterHandler &prm, const int dimen)
:
mpi_communicator (MPI_COMM_WORLD),
triangulation (mpi_communicator,
               typename Triangulation<dim>::MeshSmoothing
               (Triangulation<dim>::smoothing_on_refinement |
                Triangulation<dim>::smoothing_on_coarsening)),
dof_handler (triangulation),
ncells(prm.get_integer("number of cells in one direction")),
n_azi(prm.get_integer("quadrature order")),
ngroup(prm.get_integer("number of groups")),
is_eigen_problem(prm.get_integer("wheter (1) or not (0) to do eigenvalue calculations")),
do_nda(prm.get_integer("wheter (1) or not (0) to do NDA")),
p_order(prm.get_integer("polynomial degree")),
pcout (std::cout,
       (Utilities::MPI::this_mpi_process(mpi_communicator)
        == 0)),
computing_timer (mpi_communicator,
                 pcout,
                 TimerOutput::summary,
                 TimerOutput::wall_times)
{}

template <int dim>
EP_SN<dim>::~EP_SN()
{
  dof_handler.clear();
  delete fe;
}

template <int dim>
static void EP_SN<dim>::declare_material_parameters(ParameterHandler &prm, ParameterHandler &prm2)
{
  int nmat = prm.get_integer ("number of material types");
  int ngrp = prm.get_integer ("number of groups");
  int is_eigen = prm.get_integer ("wheter (1) or not (0) to do eigenvalue calculations");
  std::vector<unsigned int> fissile_mats;
    
  bool read_mesh = prm.get_integer ("Read in mesh or not");
  
  prm2.enter_subsection("Multigroup sigma_t, g=1 to G")
  {
    for (unsigned int m=0; m<nmat; ++m)
    {
      std::ostringstream os << "material " << m + 1;
      prm2.declare_entry(os.str(), "", Patterns::List(Patterns::Double()), "");
    }
  }
  prm2.leave_subsection();
  
  if (is_eigen)
  {
    std::ostringstream os_fiss = "list for whether materials are fissiles (0 or 1)";
    std::vector<std::string> fiss_strings = Utilities::split_string_list (prm2.get (os_fiss.str ()));
    std::ostringstream fiss_err1 << "Make sure input " << nmat << " entries for whether mats are fissile";
    AssertThrow (fiss_strings.size () == nmat,
                 ExcMessage (fiss_err1.str ()));
    for (unsigned int m=1; m<nmat; ++m)
      fissile_mats.push_back (std::atoi (fiss_strings[m].c_str()));
    std::vector<unsigned int>::iterator it = std::max_element (fissile_mats);
    AssertThrow (*it == 1,
                 ExcMessage ("No fissile material presents"));

    prm2.enter_subsection ("ksi, g=1 to G");
    {
      for (unsigned int m=0; m<nmat; ++m)
      {
        std::ostringstream os << "material " << m + 1;
        prm2.declare_entry(os.str(), "", Patterns::List(Patterns::Double()), "");
      } 
    }

    prm2.enter_subsection ("Multigroup nu*sigma_f, g=1 to G");
    {
      for (unsigned int m=0; m<nmat; ++m)
      {
        std::ostringstream os << "material " << m + 1;
        prm2.declare_entry(os.str(), "", Patterns::List(Patterns::Double()), "");
      } 
    }
    prm2.leave_subsection ();
  }

  for (unsigned int m=0; m<nmat; ++m)
  {
    std::ostringstream os << "Transfer matrix for material " << m + 1;
    prm2.enter_subsection(os.str());
    {
      for (unsigned int gin=0; gin<ngrp; ++gin)
      {
        std::ostringstream osg << "g_in " << gin + 1;
        prm2.declare_entry(osg.str(), "", Patterns::List(Patterns::Double()), "");
      }
    }
    prm2.leave_subsection();
  }

  if (is_eigen)
  {
    for (unsigned int m=0; m<nmat; ++m)
    {
      std::ostringstream osm << "Transfer matrix for material " << m + 1;
      prm2.enter_subsection(osm.str())
      {
        for (unsigned int gin=0; gin<ngrp; ++gin)
        {
          std::ostringstream osg << "g_in " << gin + 1;
          prm2.declare_entry(osg.str(), "", Patterns::List(Patterns::Double()), "");
        }
      }
      prm2.leave_subsection();
    }
  }
  
  if (!read_mesh)
  {
    prm2.subsection("Material ID map")
    {
      prm2.declare_entry("assembly material ids", "", Patterns::List(Patterns::Integer()), "Give material IDs for all blocks");
    }
    prm2.leave_subsection();
  }
}

template <int dim>
static void EP_SN<dim>::declare_parameters(ParameterHandler &prm)
{
  prm.declare_entry("number of cells in one direction", "1", Patterns::Integer(), "Number of initial cells, must be larger than 0.");
  prm.declare_entry("quadrature order", "4", Patterns::Integer(), "Gauss-Chebyshev level-symmetric-like quadrature");
  prm.declare_entry("number of groups", "1", Patterns::Integer(), "Number of groups in MG calculations");
  prm.declare_entry("wheter (1) or not (0) to do eigenvalue calculations", "1", Patterns::Integer(), "Boolean to determine problem type");
  prm.declare_entry ("list for whether materials are fissiles (0 or 1)", "", Patterns::List(Patterns::Integer()), "Whether materials are fissile");
  prm.declare_entry("wheter (1) or not (0) to do NDA", "1", Patterns::Integer(), "Boolean to determine NDA or not");
  prm.declare_entry("boundary condition types", "", Patterns::List(Patterns::Integer()), "boundary conditions types: reflective (1) and natural (0)");
  prm.declare_entry("polynomial degree", "1", Patterns::Integer(), "polynomial degree p for finite element");
  prm.declare_entry("x, y, z max values of boundary locations", "", Patterns::List(Patterns::Double()), "xmax, ymax, zmax of the boundaries, mins are zero");
  prm.declare_entry("number of cells for x, y, z directions", "", Patterns::List(Patterns::Integer()), "Geotry is hyper rectangle defined by how many cells exist per direction");
  prm.declare_entry("if boundaries are reflective", "", Patterns::List(Patterns::Integer()), "a list of 0/1 integers for all boundaries to determine if they are reflective");
  prm.declare_entry("number of material types", "", Patterns::Integer(), "must be a positive integer");
  prm.declare_entry("read in mesh or not", "0", Patterns::Integer(), "if 0, generate block-based Cartesian mesh");
  prm.declare_entry("use explicit reflective boundary condition or not", "1", Patterns::Integer(), "");
}

template <int dim>
void EP_SN<dim>::process_input_xsec ()
{
  // This block takes in sigts
  prm2.enter_subsection ("Multigroup sigma_t, g=1 to G");
  {
    for (unsigned int m=0; m<nmat; ++m)
    {
      std::ostringstream os_sigt << "material " << m + 1;
      std::vector<std::string> sigt_strings = Utilities::split_string_list (prm2.get (os_sigt.str()));
      AssertThrow (sigt_strings.size() == ngroup,
                   ExcMessage ("Ngroup is not equal to group number of sigma_t"));
      std::vector<double> tmp_sigt;
      for (unsigned int g=0; g<ngroup; ++g)
        tmp_sigt.push_back (std::atof (tmp_sigt[i].c_str ()));
      all_sigts.push_back (tmp_sigt);
    }
  }
  prm2.leave_subsection ();

  // This block takes in scattering transfer matrices
  for (unsigned int m=0; m<nmat; ++m)
  {
    std::ostringstream osm << "Transfer matrix for material " << m + 1;
    prm2.enter_subsection (osm.str ());
    {
      std::vector<std::string> sigs_strings = Utilities::split_string_list (prm2.get (os_sigs.str ()));
      std::ostringstream os_err << "Make sure input Ngroup X Ngroup entries for scattering transfer matrix for Material " << m + 1;
      AssertThrow (sigs_strings.size () == ngroup * ngroup,
                   ExcMessage (os_err.str ()));
      Table<2, double> tmp_sigs (ngroup, ngroup);
      for (unsigned int gin=0; gin<ngroup; ++gin)
      {
        for (unsigned int g=0; g<ngroup; ++g)
        {
          tmp_sigs[gin][g] = std::atof (tmp_sigs[g].c_str ());
        }
      }
    }
    prm2.leave_subsection ();
  }

  // This block is for fission
  if (is_eigen_problem)
  {
    
  }
  
}

template <int dim>
void get_cell_sigts (unsigned int &material_id, std::vector<double> &local_sigts)
{
  local_sigts = all_sigts[material_id];
}

template <int dim>
void get_cell_mfps (typename Triangulation<dim>::cell_iterator &cell,
                    std::vector<double> &local_sigts, std::vector<double> &local_mfps)
{
  double h = cell->diameter () / std::sqrt (2.0);
  for (unsigned int i=0; i<local_sigts.size (); ++i)
    local_mfps[i] = local_sigts[i] * h;
}

template <int dim>
void EP_SN<dim>::initialize_component_index ()
{
  unsigned int ind = 0;
  for (unsigned int g=0; g<ngroup; ++g)
    for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
    {
      ind += 1;
      std::pair<unsigned int, unsigned int> key (i_dir, g);
      component_index.insert (std::make_pair (key, ind));
    }
}

template <int dim>
unsigned int EP_SN<dim>::get_component_index (unsigned int &incident_angle_index, 
                                              unsigned int &g)
{
  return component_index[std::make_pair (incident_angle_index, g)];
}

template <int dim>
unsigned int EP_SN<dim>::get_reflective_direction_index (unsigned int &boundary_id, 
                                                         unsigned int &incident_angle_index)
{
  return reflective_direction_index[std::make_pair (boundary_id, incident_angle_index)];
}

void EP_SN<dim>::initialize_component_index ()
{
  unsigned int ind = 0;
  for (unsigned int g=0; g<ngroup; ++g)
    for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
    {
      ind += 1;
      std::pair<unsigned int, unsigned int> key (i_dir, g);
      component_index.insert (std::make_pair (key, ind));
    }
}

template <int dim>
unsigned int EP_SN<dim>::get_component_index (unsigned int &incident_angle_index, 
                                              unsigned int &g)
{
  return component_index[std::make_pair (incident_angle_index, g)];
}

template <int dim>
unsigned int EP_SN<dim>::get_reflective_direction_index (unsigned int &boundary_id, 
                                                         unsigned int &incident_angle_index)
{
  return reflective_direction_index[std::make_pair (boundary_id, incident_angle_index)];
}

template <int dim>
void EP_SN<dim>::initialize_ref_bc_index ()
{
  std::vector<Tensor<1, dim> > boundary_normal_vectors;
  if (dim==2)
  {
    boundary_normal_vectors.resize(4);
    boundary_normal_vectors[0][0] = -1.0;
    boundary_normal_vectors[1][0] = 1.0;
    boundary_normal_vectors[2][1] = -1.0;
    boundary_normal_vectors[3][1] = 1.0;
  }
  if (dim==3)
  {
    boundary_normal_vectors.resize(6);
    boundary_normal_vectors[0][0] = -1.0;
    boundary_normal_vectors[1][0] = 1.0;
    boundary_normal_vectors[2][1] = -1.0;
    boundary_normal_vectors[3][1] = 1.0;
    boundary_normal_vectors[4][2] = -1.0;
    boundary_normal_vectors[5][2] = 1.0;
  }
  for (unsigned int i=0; i<2*dim; ++i)
    for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
    { 
      Tensor<1, dim> out_angle = (omega_i[i_dir]
                                  *
                                  (1.0 - 2.0 * (boundary_normal_vectors[i] * omega_i[i_dir])));
      for (unsigned int r_dir=0; r_dir<n_dir; ++r_dir)
      {
        Tensor<1, dim> d_dir = out_angle;
        Tensor<1, dim> d_minus_dir = out_angle;
        d_dir -= omega_i[r_dir];
        d_minus_dir += omega_i[r_dir];
        if (d_dir.norm ()<1.0e-13 || d_minus_dir.norm ()<1.0e-13)
          reflective_direction_index.insert (std::make_pair (std::make_pair (i, i_dir), r_dir));
      }
    }
}

template <int dim>
void EP_SN<dim>::setup_system()
{
  TimerOutput::Scope t(computing_timer, "setup HO system");
  
  fe = new FE_DGQ<dim>(p_order);
  
  dof_handler.distribute_dofs (*fe);
  
  local_dofs = dof_handler.locally_owned_dofs();
  
  {
    DynamicSparsityPattern dsp(local_dofs);
    
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    
    // be careful with the following
    SparsityTools::distribute_sparsity_pattern(dsp,
                                               dof_handler.n_locally_owned_dofs_per_processor(),
                                               mpi_communicator,
                                               local_dofs);
  }
  
  for (unsigned int g=0; g<ngroup; ++g)
  {
    if (do_nda)
    {
      vec_lo_sys.push_back(new LA::MPI::SparseMatrix);
      vec_lo_rhs.push_back(new LA::MPI::Vector);
      vec_lo_sflx.push_back(new LA::MPI::Vector);
      vec_lo_sflx_old.push_back(new LA::MPI::Vector);
    }
    
    vec_ho_rhs.push_back(new LA::MPI::Vector);
    vec_ho_sflx.push_back(new LA::MPI::Vector);
    vec_ho_sflx_old.push_back(new LA::MPI::Vector);
    
    for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
    {
      vec_ho_sys.push_back(new LA::MPI::SparseMatrix);
      vec_aflx.push_back(new LA::MPI::Vector);
    }
  }
  
  for (unsigned int g=0; g<ngroup; ++g)
  {
    if (do_nda)
    {
      vec_lo_sys[g]->reinit(local_dofs,
                            local_dofs,
                            dsp,
                            mpi_communicator);
      vec_lo_rhs[g]->reinit(local_dofs,
                            mpi_communicator);
      vec_lo_sflx[g]->reinit(local_dofs,
                             mpi_communicator);
      vec_lo_sflx_old[g]->reinit(local_dofs,
                                 mpi_communicator);
    }
    
    vec_ho_rhs[g]->reinit(local_dofs,
                          mpi_communicator);
    vec_ho_sflx[g]->reinit(local_dofs,
                           mpi_communicator);
    vec_ho_sflx_old[g]->reinit(local_dofs,
                               mpi_communicator);
    
    for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
    {
      vec_ho_sys[get_component_index(i_dir, g)]->reinit(local_dofs,
                                                        local_dofs,
                                                        dsp,
                                                        mpi_communicator);
      vec_aflx[get_component_index(i_dir, g)]->reinit(local_dofs,
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
 
 for (unsigned int g=0; g<ngroup; ++g)
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
void EP_SN<dim>::setup_boundary_types(ParameterHandler &prm)
{
  std::vector<std::string> boundary_locations = Utilities::split_string_list (
      prm.get("x, y, z max and min values of boundary locations"));
  AssertThrow (boundary_locations.size()==2*dim, 
      ExcMessage("List of boundary location specifications should be identical to twice the dimension"))
  for (unsigned int i=0; i<boundary_locations.size(); ++i)
    boundary_location_values.push_back(std::atof(boundary_locations[i].c_str()));
  
  typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active (), endc = triangulation.end ();
  for (; cell!=endc; ++cell)
  {
    if (cell->is_locally_owned)
    {
      for (unsigned int fn=0; fn<GeometryInfo<dim>::face_per_cell; ++fn)
      {
        if (cell->face(fn)->at_boundary())
        {
          Point<dim> ct = cell->face(fn)->center();
          // left boundary
          if (std::fabs(ct[0]-xmin)<1.0e-14)
            cell->face(fn)->set_boundary_id(0);
          
          // right boundary
          if (std::fabs(ct[0]-xmax)<1.0e-14)
            cell->face(fn)->set_boundary_id(1);
          
          // 2D and 3D boundaries
          if (ct.size()>1)
          {
            // 2D boundaries
            // front boundary
            if (std::fabs(ct[1]-ymin)<1.0e-14)
              cell->face(fn)->set_boundary_id(2);
            
            // rear boundary
            if (std::fabs(ct[1]-ymax)<1.0e-14)
              cell->face(fn)->set_boundary_id(3);
            
            // 3D boundaries
            if (ct.size()>2)
            {
              // front boundary
              if (std::fabs(ct[2]-zmin)<1.0e-14)
                cell->face(fn)->set_boundary_id(4);
              
              // rear boundary
              if (std::fabs(ct[2]-zmax)<1.0e-14)
                cell->face(fn)->set_boundary_id(5);
            }
          }
        }
      }// face
    }// locally owned cell
  }// cell
  
  std::vector<std::string> ref_bc_bools = Utilities::split_string_list(prm.get("if boundaries are reflective"));
  AssertThrow(ref_bc_bools.size()==2*dim, ExcMessage("List of boundary type specifications should be identical to twice the dimension"))
  for (unsigned int i=0; i<ref_bc_bools.size(); ++i)
    is_reflective_bc.push_back(Utilities::string_to_int(ref_bc_bools[i].c_str()));
}

template <int dim>
void EP_SN<dim>::assemble_ho_system ()
{
  TimerOutput::Scope t(computing_timer, "assembly HO");
  
  const QGauss<dim>  q_rule(p_order+1);
  const QGauss<dim>  qf_rule(p_order+1);
  
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
  // face finite element object for the side of the face in neighbor cell
  FEFaceValues<dim> fvf_nei(*fe, qf_rule,
                            update_values | update_gradients |
                            update_quadrature_points | update_normal_vectors |
                            update_JxW_values);
  
  
  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q = q_rule.size();
  const unsigned int n_qf = qf_rule.size();
  
  // volumetric matrix
  std::vector<FullMatrix<double> > cell_matrix(n_dir * ngroup, FullMatrix<double>(dofs_per_cell, dofs_per_cell));
  FullMatrix<double> cell_matrix_collision(dofs_per_cell, dofs_per_cell);
  
  // face terms: v^\pm * u^\pm
  std::vector<FullMatrix<double> > all_real_vp_up(ngroup * n_dir, FullMatrix<double>(dofs_per_cell, dofs_per_cell));
  std::vector<FullMatrix<double> > all_real_vp_un(ngroup * n_dir, FullMatrix<double>(dofs_per_cell, dofs_per_cell));
  std::vector<FullMatrix<double> > all_real_vn_up(ngroup * n_dir, FullMatrix<double>(dofs_per_cell, dofs_per_cell));
  std::vector<FullMatrix<double> > all_real_vn_un(ngroup * n_dir, FullMatrix<double>(dofs_per_cell, dofs_per_cell));
  
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  std::vector<types::global_dof_index> neigh_dof_indices (dofs_per_cell);
  
  // volumetric pre-assembly matrices
  std::vector<FullMatrix<double> > vec_collision_local (n_q, FullMatrix<double>(dofs_per_cell, dofs_per_cell));
  std::vector<std::vector<FullMatrix<double> > > vec_streaming_local (n_q, 
      std::vector<FullMatrix<double> > (n_dir, FullMatrix<double> (dofs_per_cell, dofs_per_cell)));
  
  // face pre-assembly matrices: value penalty
  std::vector<std::vector<FullMatrix<double> > > vec_vp_up (n_q, 
      FullMatrix<double> (dofs_per_cell, dofs_per_cell));
  std::vector<std::vector<FullMatrix<double> > > vec_vp_un (n_q, 
      FullMatrix<double> (dofs_per_cell, dofs_per_cell));
  std::vector<std::vector<FullMatrix<double> > > vec_vn_up (n_q, 
      FullMatrix<double> (dofs_per_cell, dofs_per_cell));
  std::vector<std::vector<FullMatrix<double> > > vec_vn_un (n_q, 
      FullMatrix<double> (dofs_per_cell, dofs_per_cell));
  
  // face pre-assembly matrices: gradient penalty 1
  std::vector<std::vector<FullMatrix<double> > > vec_dvp_up (n_q, 
      std::vector<FullMatrix<double> > (n_dir, FullMatrix<double> (dofs_per_cell, dofs_per_cell)));
  std::vector<std::vector<FullMatrix<double> > > vec_dvp_un (n_q, 
      std::vector<FullMatrix<double> > (n_dir, FullMatrix<double> (dofs_per_cell, dofs_per_cell)));
  std::vector<std::vector<FullMatrix<double> > > vec_dvn_up (n_q, 
      std::vector<FullMatrix<double> > (n_dir, FullMatrix<double> (dofs_per_cell, dofs_per_cell)));
  std::vector<std::vector<FullMatrix<double> > > vec_dvn_un (n_q, 
      std::vector<FullMatrix<double> > (n_dir, FullMatrix<double> (dofs_per_cell, dofs_per_cell)));
  
  // face pre-assembly matrices: gradient penalty 2
  std::vector<std::vector<FullMatrix<double> > > vec_vp_dup (n_q, 
      std::vector<FullMatrix<double> > (n_dir, FullMatrix<double> (dofs_per_cell, dofs_per_cell)));
  std::vector<std::vector<FullMatrix<double> > > vec_vp_dun (n_q, 
      std::vector<FullMatrix<double> > (n_dir, FullMatrix<double> (dofs_per_cell, dofs_per_cell)));
  std::vector<std::vector<FullMatrix<double> > > vec_vn_dup (n_q, 
      std::vector<FullMatrix<double> > (n_dir, FullMatrix<double> (dofs_per_cell, dofs_per_cell)));
  std::vector<std::vector<FullMatrix<double> > > vec_vn_dun (n_q, 
      std::vector<FullMatrix<double> > (n_dir, FullMatrix<double> (dofs_per_cell, dofs_per_cell)));
  
  int pre_assemble_cell = 0;
  int pre_assemble_face = 0;
  
  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (typename DoFHandler<dim>::active_cell_iterator cell=dof_handler.begin_active(); 
      cell!=dof_handler.end(); ++cell)
  {
    if (cell->is_locally_owned())
    {
      fv.reinit(cell);
      cell->get_dof_indices(local_dof_indices);
      std::vector<double> local_sigts (ngroup);
      std::vector<double> local_mfps (ngroup);
      get_cell_sigts (cell, local_sigts);
      get_cell_mfps (cell, local_sigts, local_mfps);
      // FixIt: a more proper definition for h_cell
      double cell_mfp = cell_sigt * cell->diameter() / std::sqrt(2.0);
      
      if (pre_assemble_cell==0)
      {
        for (unsigned int qi=0; qi<n_q; ++qi)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              vec_collision_local[qi](i,j) += (fv.shape_value(i,qi) *
                                               fv.shape_value(j,qi));
        
        
        for (unsigned int qi=0; qi<n_q; ++qi)
          for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              for (unsigned int j=0; j<dofs_per_cell; ++j)
                vec_streaming_local[qi][i_dir](i,j) += ((fv.shape_grad(i,qi) *
                                                         omega_i[i_dir]) *
                                                        (fv.shape_grad(j,qi) *
                                                         omega_i[i_dir]));
        
        pre_assemble_cell = 1;
      }
      
      // Using mass_matrix would be benificial for multigroup assembly
      FullMatrix<double> mass_matrix(dofs_per_cell, dofs_per_cell);
      std::vector<FullMatrix<double> > stiffness(n_dir, FullMatrix<double>(dofs_per_cell, dofs_per_cell));
      
      for (unsigned int g=0; g<ngroup; ++g)
      {
        if (g==0)
          for (unsigned int qi=0; qi<n_q; ++qi)
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              for (unsigned int j=0; j<dofs_per_cell; ++i)
                mass_matrix(i,j) += vec_collision_local[qi](i,j) * fv.JxW(qi);
        
        // cell_matrix_collision = mass_matrix
        // cell_matrix_collision *= cell_sigt;
        
        for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
        {
          if (g==0)
            for (unsigned int qi=0; qi<n_q; ++qi)
              for (unsigned int i=0; i<dofs_per_cell; ++i)
                for (unsigned int j=0; j<dofs_per_cell; ++j)
                  stiffness[i_dir](i,j) += vec_streaming_local[qi][i_dir](i,j) * fv.JxW(qi);
          
          int ind = index(i_dir, g);
          cell_matrix[ind] = stiffness[i_dir];
          cell_matrix[ind] /= cell_sigt;
          // cell_matrix[ind].add(1.0, cell_matrix_collision);
          cell_matrix[ind].add(cell_sigt, mass_matrix);
          
          /*
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              vec_ho_sys[ind](local_dof_indices[i],local_dof_indices[j]) += cell_matrix(i,j);
           */
        }// i_dir
      }// g
      
      for (unsigned int fn=0; fn<GeometryInfo<dim>::face_per_cell; ++fn)
      {
        // BC imposition: only vacuum BCs need to be imposed in bilinear form
        // Reflective boundary will only exist in linear form in source iterations
        if (cell->face(fn)->at_boundary())
        {
          fvf.reinit(cell, fn);
          Tensor<1, dim> vec_n = fvf.get_normal_vectors()[0];
          unsigned int boundary_id = cell->face(fn)->boundary_id ();
          if (!is_reflective_bc[boundary_id])
          {
            for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
            {
              double absndo;
              // Note: we assume the face is not curvilinear
              absndo = std::fabs(vec_n * omega_i[i_dir]);
              for (unsigned int g=0; g<ngroup; ++g)
              {
                int ind = index(i_dir, g);
                for (unsigned int qi=0; qi<n_qf; ++qi)
                  for (unsigned int i=0; i<dofs_per_cell; ++i)
                    for (unsigned int j=0; j<dofs_per_cell; ++j)
                      cell_matrix[ind](i,j) += (absndo *
                                                fvf.shape_value(i,qi) *
                                                fvf.shape_value(j,qi) *
                                                fvf.JxW(qi));
              }// g
            }// i_dir
          }
          else
          {
            if (is_explicit_reflective)// assemble nothing if false
              for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
              {
                unsigned int r_dir = get_reflective_direction_index (boundary_id, i_dir);
                double absndo = omega_i[i_dir] * vec_n;
                for (unsigned int g=0; g<ngroup; ++g)
                {
                  unsigned int ind = get_component_index (i_dir, g);
                  for (unsigned int qi=0; qi<n_qf; ++qi)
                    for (unsigned int i=0; i<dofs_per_cell; ++i)
                      for (unsigned int j=0; j<dofs_per_cell; ++j)
                        cell_matrix[ind](i,j) += -(absndo *
                                                   fvf.shape_value(i,qi) *
                                                   (omega_i[r_dir] * fvf.shape_grad(j,qi)) /
                                                   local_sigts[g] *
                                                   fvf.JxW(qi));
                }// g
              }// i_dir
          }// is_reflective_bc
        }// boundary faces for robin boundaries
        
        if (!cell->face(fn)->at_boundary() &&
            cell->id()>cell->neighbor(fn)->id())
        {
          typename DoFHandler<dim>::cell_iterator neigh = cell->neighbor(fn);
          // initialize the elements of sides of the face in current cell and neighbor
          fvf.reinit(cell, fn);
          fvf_nei.reinit(neigh, cell->neighbor_face_no(fn));
          Tensor<1,dim> n_vec = fvf.get_normal_vectors()[0];
          double sige;
          std::vector<double> neigh_sigts;
          std::vector<double> neigh_mfps;
          get_cell_sigts (neigh, neigh_sigts);
          get_cell_mfps (neigh, neigh_sigts, neigh_mfps);
          
          if (pre_assemble_face==0)
          {
            // assemble once for vp/n, up/n in a reference cell
            for (unsigned int qi=0; qi<n_qf; ++qi)
            {
              double jxw = fvf.JxW(qi);
              for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
                for (unsigned int g=0; g<ngroup; ++g)
                {
                  int ind = index(i_dir, g);
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
            
            pre_assemble_face = 1;
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
                  for (unsigned int g=0; g<ngroup; ++g)
                  {
                    int ind = index(i_dir, g);
                    // get cross sections ones per cell/neighbor
                    if (qi==0 && i==0 && j==0)
                    {
                      if (i_dir==0)
                      {
                        neigh_sigt = sigt[neigh->material_id()][g];
                        neigh_mfp = neigh_sigt * neigh->diameter() / std::sqrt(2.0);
                      }
                      sige = std::max(0.25, (tensor_norms[i_dir] / cell_mfp + tensor_norms[i_dir] / neigh_mfp));
                    }
                    
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
                                                 dvp_up_jxw / cell_sigt
                                                 -
                                                 vp_dup_jxw / cell_sigt);
                    
                    all_real_vp_un[ind](i,j) += (-sige * vp_un_jxw
                                                 +
                                                 dvp_un_jxw / cell_sigt
                                                 -
                                                 vp_dun_jxw / neigh_sigt);
                    
                    all_real_vn_up[ind](i,j) += (-sige * vp_up_jxw
                                                 -
                                                 dvn_up_jxw / neigh_sigt
                                                 +
                                                 vn_dup_jxw / cell_sigt);
                    
                    all_real_vn_un[ind](i,j) += (sige * vp_up_jxw
                                                 +
                                                 dvn_un_jxw / neigh_sigt
                                                 +
                                                 vn_dun_jxw / neigh_sigt);
                  }// g
                }// i_dir
              }// j
          }// qi
          
          neigh->get_dof_indices(neigh_dof_indices);
          
          for (unsigned int k=0; k<n_total_ho_vars; ++k)
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              for (unsigned int j=0; j<dofs_per_cell; ++j)
              {
                vec_ho_sys[k](local_dof_indices[i],local_dof_indices[j]) += all_real_vp_up[k](i,j);
                
                vec_ho_sys[k](local_dof_indices[i],neigh_dof_indices[j]) += all_real_vp_un[k](i,j);
                
                vec_ho_sys[k](neigh_dof_indices[i],local_dof_indices[j]) += all_real_vn_up[k](i,j);
                
                vec_ho_sys[k](neigh_dof_indices[i],neigh_dof_indices[j]) += all_real_vn_un[k](i,j);
              }// j
        }// non-boundary face
      }// face
      
      for (unsigned int k=0; k<n_total_ho_vars; ++k)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            vec_ho_sys[k](local_dof_indices[i],local_dof_indices[j]) += cell_matrix[k](i,j);
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
void EP_SN<dim>::angular_quad ()
{
  Assert (n_azi%2==0, ExcDimensionMismatch(n_azi%2, 1));
  
  QGauss<1> mu_quad (n_azi);
  
  std::ofstream quadr;
  quadr.open("quadr.txt");
  quadr << "Dim = " << dim << ", SN order = " << n_azi << std::endl;
  
  if (dim==1)
  {
    total_angle = 2.0;
    quadr << "1D mu and weights:" << std::endl;
    for (unsigned int i=0; i<n_azi/2; ++i)
    {
      Tensor <1, dim> tmp;
      tmp[0] = mu_quad.point(i)[0] * 2.0 - 1.0;;
      omega_i.push_back (tmp);
      wi.push_back (mu_quad.weight(i))
      quadr << tmp[0] << "," << wi[i] << std::endl;
    }
    n_dir = n_azi / 2;
  }
  
  if (dim==2)
  {
    total_angle = 4.0 * pi;
    quadr << "2D Omega_i and weights:" << std::endl;
    
    for (unsigned int i=n_azi/2; i<n_azi; ++i)
    {
      Tensor<1, dim> tmp;
      double mut = mu_quad.point(i)[0]*2.0 - 1.0;
      int level_angle_num = 4 * (n_azi - i);
      double delta_level_angle = 2.0 * pi / level_angle_num;
      double level_weight = mu_quad.weight(i) * 4.0 * pi;
      
      for (unsigned int j=0; j<level_angle_num; ++j)
      {
        double angle = 0.5 * delta_level_angle + (double)(j) * delta_level_angle;
        tmp[0] = std::sqrt (1.0 - mut * mut) * cos (angle);
        tmp[1] = std::sqrt (1.0 - mut * mut) * sin (angle);
        // solely for EP quadratures in 2D.
        if (tmp[0]>0)
        {
          omega_i.push_back(tmp);
          double point_wt = level_weight / level_angle_num;
          wi.push_back(point_wt);
          quadr << tmp[0] << ", " << tmp[1] << ", " << mut << ", " << point_wt << std::endl;
        }
        
      }
      
    }
    n_dir = n_azi * (n_azi + 2) / 4;
  }
  
  if (dim==3)
  {
    total_angle = 4.0 * pi;
    quadr << "3D Omega_i and weights:" << std::endl;
    
    for (unsigned int i=0; i<n_azi/2; ++i)
    {
      int level_angle_num = i < n_azi / 2 ? 4 * (i + 1) : 4 * (n_azi - i);
      double delta_level_angle = 2.0 * pi / level_angle_num;
      double level_weight = mu_quad.weight(i) * 4.0 * pi;
      Tensor<1, dim> tmp;
      double mut = mu_quad.point(i)[0] * 2.0 - 1.0;
      
      for (unsigned int j=0; j<level_angle_num; ++j)
      {
        double angle = 0.5 * delta_level_angle + (double)(j) * delta_level_angle;
        tmp[0] = std::sqrt(1.0 - mut * mut) * cos(angle);
        tmp[1] = std::sqrt(1.0 - mut * mut) * sin(angle);
        tmp[2] = mut;
        omega_i.push_back(tmp);
        double point_wt = level_weight / level_angle_num;
        wi.push_back(point_wt);
        quadr << tmp[0] << ", " << tmp[1] << ", " << tmp[2] << ", " << point_wt << std::endl;
      }
      
    }
    n_dir = n_azi * (n_azi + 2) / 4;
  }
  n_total_ho_vars = n_dir * ngroup;
  quadr.close();
  
  for (unsigned int i=0; i<n_total_ho_vars; ++i)
  {
    FEValuesExtractors::Scalar tmp(i);
    comp.push_back(tmp);
  }
  
  // estimate tensor norm to do penalty method
  for (unsigned int i=0; i<n_dir; ++i)
  {
    Tensor<2, dim> tensor_tmp = outer_product(omega_i[i_dir], omega_i[i_dir]);
    tensor_norms.push_back(tensor_tmp.norm());
  }
}

template <int dim>
void EP_SN<dim>::ho_solve (unsigned int &n_iter)
{
  TimerOutput::Scope t(computing_timer, "HO solve");
  LA::MPI::Vector
  completely_distributed_solution (ho_local_dofs,
                                   mpi_communicator);
  
  SolverControl solver_control (ho_dof_handler.n_dofs(), ho_rhs.l2_norm() * 1e-12);
  
#ifdef USE_PETSC_LA
  LA::SolverCG solver(solver_control, mpi_communicator);
#else
  LA::SolverCG solver(solver_control);
#endif
  
  // LA::MPI::PreconditionAMG preconditioner;
  
  if (ho_precond_kind == "amg")
  {
    if (n_iter==0)
    {
      pre_AMG.reset();
      pre_AMG = (std_cxx11::shared_ptr<LA::MPI::PreconditionAMG> (new LA::MPI::PreconditionAMG));
      LA::MPI::PreconditionAMG::AdditionalData data;
      
#ifdef USE_PETSC_LA
      data.symmetric_operator = true;
#else
      /* Trilinos defaults are good */
#endif
      pre_AMG->initialize(ho_mat, data);
    }
    
    //    solver.solve (system_matrix, completely_distributed_solution, system_rhs, preconditioner);
    solver.solve (ho_sys, all_aflx, ho_rhs, pre_AMG);
  }
  // FixIt: add other preconditioners
  // constraints.distribute (completely_distributed_solution);
}

template <int dim>
void EP_SN<dim>::ho_solve_bicgstab (unsigned int &n_iter)
{
  TimerOutput::Scope t(computing_timer, "HO solve");
  LA::MPI::Vector
  completely_distributed_solution (ho_local_dofs,
                                   mpi_communicator);
  
  SolverControl solver_control (ho_dof_handler.n_dofs(), ho_rhs.l2_norm() * 1e-12);
  
#ifdef USE_PETSC_LA
  LA::SolverBicgstab solver(solver_control, mpi_communicator);
#else
  LA::SolverBicgstab solver(solver_control);
#endif
  
  // LA::MPI::PreconditionAMG preconditioner;
  
  if (ho_precond_kind == "amg")
  {
    if (n_iter==0)
    {
      pre_AMG.reset();
      pre_AMG = (std_cxx11::shared_ptr<LA::MPI::PreconditionAMG> (new LA::MPI::PreconditionAMG));
      LA::MPI::PreconditionAMG::AdditionalData data;
      
#ifdef USE_PETSC_LA
      data.symmetric_operator = false;
#else
      /* Trilinos defaults are good */
#endif
      pre_AMG->initialize(ho_mat, data);
    }
    
    //    solver.solve (system_matrix, completely_distributed_solution, system_rhs, preconditioner);
    solver.solve (ho_sys, all_aflx, ho_rhs, pre_AMG);
  }
  // FixIt: add other preconditioners
  // constraints.distribute (completely_distributed_solution);
}

/*
template <int dim>
void EP_SN<dim>::lo_solve()
{
  TimerOutput::Scope t(computing_timer, "LO solve");
  LA::MPI::Vector
  completely_distributed_solution (ho_local_dofs,
                                   mpi_communicator);
  
  SolverControl solver_control (ho_dof_handler.n_dofs(), ho_rhs.l2_norm() * 1e-12);
  
#ifdef USE_PETSC_LA
  LA::SolverCG solver(solver_control, mpi_communicator);
#else
  LA::SolverCG solver(solver_control);
#endif
  
  // LA::MPI::PreconditionAMG preconditioner;
  
  if (lo_precond_kind == "amg")
  {
    for (unsigned int g=0; g<ngroup; ++g)
    {
      pre_AMG.reset();
      pre_AMG = (std_cxx11::shared_ptr<LA::MPI::PreconditionAMG> (new LA::MPI::PreconditionAMG));
      LA::MPI::PreconditionAMG::AdditionalData data;
      
#ifdef USE_PETSC_LA
      data.symmetric_operator = true;
#else

#endif
      pre_AMG->initialize(*(vec_lo_sys)[g], data);
      
      solver.solve (*(vec_lo_sys)[g], *(vec_lo_solu)[g], *(vec_lo_rhs)[g], pre_AMG);
    }
  }
  // FixIt: add other preconditioners
  // constraints.distribute (completely_distributed_solution);
}
*/

template <int dim>
void EP_SN<dim>::generate_moments ()
{
  // FitIt: only scalar flux is generated for now
  AssertThrow(do_nda==false, ExcMessage("Moments are generated only without NDA"));
  if (!do_nda)
  {
    for (unsigned int g=0; g<ngroup; ++g)
    {
      vec_ho_sflx[g] = 0;
      for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
        vec_ho_sflx[g].add(4.0*wi[i_dir], vec_aflx[index(i_dir, g)]);
    }
  }
}

template <int dim>
void EP_SN<dim>::generate_source ()
{
  const QGauss<dim>  q_rule(p_order+1);
  const QGauss<dim>  qf_rule(p_order+1);
  
  unsigned int n_q = q_rule.size();
  unsigned int n_qf = q_rule.size();
  
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
  std::vector<Vector<double> > vec_cell_rhs(ngroup, Vector<double> (dofs_per_cell));
  std::vector<Vector<double> > vec_cell_rhs_reflective_bc(n_total_ho_vars, Vector<double> (dofs_per_cell)); 
  
  for (typename DoFHandler<dim>::active_cell_iterator cell = triangulation.begin_active(); cell!= triangulation.end(); ++cell)
  {
    if (cell->is_locally_owned())
    {
      fv.reinit(cell);
      unsigned int material_id = cell->material_id ();
      std::vector<std::vector<double> > all_cell_sflx(ngroup, std::vector<double>(n_q));
      for (unsigned int gin=0; gin<ngroup; ++gin)
      {
        if (do_nda)
        {
          // FixIt
        }
        else
          fv.get_function_values(vec_ho_sflx[gin], all_cell_sflx[gin]);
      }// gin
      
      std::vector<std::vector<double> > cell_transfers(ngroup, std::vector<double>(ngroup));
      for (unsigned int g=0; g<ngroup; ++g)
      {
        q_per_ster[g] = is_eigen_problem ? 0.0 : source_per_ster[material_id][g];

        for (unsigned int gin=0; gin<ngroup; ++gin)
          cell_transfers_per_ster[gin][g] = (sigs_per_ster[material_id][gin][g] 
                                             +
                                             is_eigen_problem ? 
                                             fiss_per_ster_over_k[material_id][gin][g] : 0.0);
      }

      for (unsigned int gin=0; gin<ngroup; ++gin)
        fv.get_function_values(vec_ho_sflx[gin], all_cell_sflx[gin]);
      
      for (unsigned int qi=0; qi<n_q; ++qi)
      {
        // do something
        double jxw = fv.JxW(qi);
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          double real_test_func = fv.shape_value(i,qi) * jxw;
          for (unsigned int g=0; g<ngroup; ++g)
            for (unsigned int gin=0; gin<ngroup; ++gin)
            {
              
              vec_cell_rhs[g](i) += (real_test_func * 
                                     (cell_transfers[gin][g] * all_cell_sflx[gin][qi]
                                      +
                                      is_eigen_problem ?
                                      0.0 : q_per_ster[g]));
            }
        }
      }// qi
      
      if (cell->at_boundary())
      {
        // Boundary parts
        for (unsigned int fn=0; fn<GeometryInfo<dim>::face_per_cell; ++fn)
        {
          if (cell->at_boundary (fn) &&
              (is_reflective_bc[boundary_id] && (!is_explicit_reflective)))
          {
            unsigned int boundary_id = cell->face(fn)->boundary_id ();
            fvf.reinit(cell, fn);
            for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
            {
              unsigned int r_dir = get_reflective_direction_index (boundary_id, i_dir);
              for (unsigned int g=0; g<ngroup; ++g)
              {
                int ind = index(i_dir, g);
                Tensor<1, dim> vec_n = fvf.get_normal_vectors()[0];
                std::vector<Tensor<1, dim> > cell_daflx(n_qf);
                fvf.get_function_gradients(vec_aflx[get_component_index (r_dir, g)], cell_daflx);
                for (unsigned int qi=0; qi<n_qf; ++qi)
                  for (unsigned int i=0; i<dofs_per_cell; ++i)
                    vec_cell_rhs[ind](i) += (fvf.shape_value(i, qi) *
                                             vec_n * omega_i[i_dir] / cell_sigt *
                                             omega_i[r_dir] * cell_daflx[qi] *
                                             fvf.JxW(qi));               
              }// g
            }// i_dir
          }// reflective boundary face
        }
      }
    }// local cells
  }
}

template <int dim>
void EP_SN<dim>::output_results(unsigned int g) const
{
  DataOut<dim> data_out;
  data_out.attach_dof_handler(lo_dof_handler);
  data_out.add_data_vector(lo_local_dofs, "u");
  
  Vector<float> subdomain(triangulation.n_active_cells());
  for (unsigned int i=0; i<subdomain.size(); ++i)
    subdomain(i) = triangulation.locally_owned_subdomain();
  data_out.add_data_vector(subdomain, "subdomain");
  
  data_out.build_patches();
  
  const std::string filename = ("solution-" +
                                Utilities::int_to_string (g, 2) +
                                "." +
                                Utilities::int_to_string
                                (triangulation.locally_owned_subdomain(), 4));
  std::ofstream output ((filename + ".vtu").c_str());
  data_out.write_vtu (output);
  
  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
  {
    std::vector<std::string> filenames;
    for (unsigned int i=0;
         i<Utilities::MPI::n_mpi_processes(mpi_communicator);
         ++i)
      filenames.push_back ("solution-" +
                           Utilities::int_to_string (g, 2) +
                           "." +
                           Utilities::int_to_string (i, 4) +
                           ".vtu");
    
    std::ofstream master_output (("solution-" +
                                  Utilities::int_to_string (cycle, 2) +
                                  ".pvtu").c_str());
    data_out.write_pvtu_record (master_output, filenames);
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
void EP_SN<dim>::power_iteration ()
{
}

template <int dim>
void EP_SN<dim>::source_iteration ()
{
}

template <int dim>
void EP_SN<dim>::generate_globally_refined_grid ()
{

  std::vector<std::string> cell_strings = Utilities::split_string_list (prm2.get ("number of cells for x, y, z directions"));
  AssertThrow (cell_strings.size() == dim,
               ExcMessage ("Entries for numbers of cells should be equal to dimension"));
  std::vector<unsigned int> cells_per_dir;
  for (unsigned int d=0; d<dim; ++d)
    cells_per_dir.push_back (std::atof (cell_strings[i].c_str ()));
  prm.get ();
  std::vector<std::vector<double> > spacings (dim);
  GridGenerator::subdivided_hyper_rectangle (triangulation, spacings,
                                             lower_left_point, 
                                             material_id_table);
  triangulation.refine_global (global_refinements);
}

template <int dim>
void EP_SN<dim>::report_system ()
{
  pcout << "Number of active cells: "
  << triangulation.n_global_active_cells()
  << std::endl
  << "Number of high-order degrees of freedom: "
  << ho_dof_handler.n_dofs()
  << std::endl;
    
  if (do_nda)
  {
    pcout << "Number of low-order degrees of freedom: "
    << lo_dof_handler.n_dofs() * ngroup
    << std::endl;
  }
}

template <int dim>
void EP_SN<dim>::do_iterations ()
{
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
void EP_SN<dim>::run()
{
 
  generate_globally_refined_grid ();

  setup_system ();

  report_system ();

  assemble_ho_system ();

  do_iterations ();

  if (Utilities::MPI::n_mpi_processes(mpi_communicator) <= 32)
  {
    TimerOutput::Scope t(computing_timer, "output");
    output_results();
  }

  computing_timer.print_summary ();
  computing_timer.reset ();
}

int main(int argc, char *argv[])
{
  try
  {
    using namespace dealii;
    using namespace Step40;
    
    int dimension;
    if (argc!=4)
    {
      std::cerr << "Call the program as ./dg-ep-proto input1 input2 dimension" << std::endl;
      return 1;
    }
    else
    {
      std::stringstream convert(argv[3]);
      convert >> dimension;
      assert (dimension == 2 || dimension == 3);
    }
    
    ParameterHandler prm;
    ParameterHandler xsec_prm;
    EP_SN<dimension>::declare_parameters(prm);
    prm.read_input(argv[1]);
    EP_SN<dimension>::declare_material_parameters(prm, xsec_prm);
    xsec_prm.read_input(argv[2]);
    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
    EP_SN<dimension> DG_EP (prm, xsec_prm);
    
    DG_EP.run();
  }
  catch (std::exception &exc)
  {
    std::cerr << std::endl << std::endl
    << "----------------------------------------------------"
    << std::endl;
    std::cerr << "Exception on processing: " << std::endl
    << exc.what() << std::endl
    << "Aborting!" << std::endl
    << "----------------------------------------------------"
    << std::endl;
    
    return 1;
  }
  catch (...)
  {
    std::cerr << std::endl << std::endl
    << "----------------------------------------------------"
    << std::endl;
    std::cerr << "Unknown exception!" << std::endl
    << "Aborting!" << std::endl
    << "----------------------------------------------------"
    << std::endl;
    return 1;
  }
  
  return 0;
}
