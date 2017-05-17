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
#include <deal.II/base/table.h>

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
  // using namespace dealii::LinearAlgebraTrilinos;
  using namespace dealii::LinearAlgebraPETSc;
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
#include <deal.II/fe/fe_poly.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/numbers.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

#include <boost/algorithm/string.hpp>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace dealii;

template <int dim>
class EP_SN
{
public:
  EP_SN (const ParameterHandler &prm);
  ~EP_SN ();
  
  void run (ParameterHandler &prm);
  
  static void declare_parameters (ParameterHandler &prm0);
  
private:
  void process_input (ParameterHandler &prm);
  void process_problem_definition (ParameterHandler &prm);
  void preprocess_reflective_bc (ParameterHandler &prm);
  void process_coordinate_information (ParameterHandler &prm);
  void process_material_properties (ParameterHandler &prm);
  void process_eigen_material_properties (ParameterHandler &prm);
  void initialize_relative_position_to_id_map (ParameterHandler &prm);
  void setup_system ();
  void generate_globally_refined_grid ();
  void initialize_material_id ();
  void report_system ();
  // void setup_lo_system();
  void setup_boundary_ids ();
  void get_cell_mfps (unsigned int &material_id, double &cell_dimension,
                      std::vector<double> &local_mfps);
  void assemble_ho_volume_boundary ();
  void assemble_ho_interface ();
  void assemble_ho_system ();
  void do_iterations ();
  void angular_quad ();
  void initialize_ref_bc_index ();
  void initialize_component_index ();

  unsigned int get_component_index (unsigned int &incident_angle_index, unsigned int &g);
  unsigned int get_direction (unsigned int &comp_ind);
  unsigned int get_component_group (unsigned int &comp_ind);
  unsigned int get_reflective_direction_index (unsigned int &boundary_id, 
                                               unsigned int &incident_angle_index);
  void get_cell_relative_position (Point<dim> &position,
                                   std::vector<unsigned int> &relative_position);
  
  void assemble_lo_system ();
  void prepare_correction_aflx ();
  
  void initialize_ho_preconditioners (); 
  void ho_solve ();
  void lo_solve ();
  void refine_grid ();
  void output_results () const;
  void generate_moments ();
  void generate_ho_source ();
  void generate_fixed_source ();
  void power_iteration ();
  void source_iteration ();
  void scale_fiss_transfer_matrices ();
  void renormalize_sflx (std::vector<LA::MPI::Vector*> &target_sflxes, double &normalization_factor);
  
  double estimate_k (double &fiss_source,
                     double &fiss_source_prev_gen,
                     double &k_prev_gen);
  double estimate_fiss_source (std::vector<LA::MPI::Vector*> &phis);
  double estimate_phi_diff (std::vector<LA::MPI::Vector*> &phis_newer,
                            std::vector<LA::MPI::Vector*> &phis_older);
  
  
  void NDA_PI ();
  void NDA_SI ();
  
  MPI_Comm mpi_communicator;
  
  parallel::distributed::Triangulation<dim> triangulation;
  
  DoFHandler<dim> dof_handler;
  // FE_DGQ<dim> *fe;
  FE_Poly<TensorProductPolynomials<dim>,dim,dim> *fe;
  
  // FixIt: involve relevant_dofs for future if refinement is necessary
  IndexSet local_dofs;
  
  const double pi;
  
  // HO system
  std::vector<LA::MPI::SparseMatrix*> vec_ho_sys;
  std::vector<LA::MPI::Vector*> vec_aflx;
  std::vector<LA::MPI::Vector*> vec_ho_rhs;
  std::vector<LA::MPI::Vector*> vec_ho_fixed_rhs;
  std::vector<LA::MPI::Vector*> vec_ho_sflx;
  std::vector<LA::MPI::Vector*> vec_ho_sflx_old;
  std::vector<LA::MPI::Vector*> vec_ho_sflx_prev_gen;
  double k_ho;
  double k_ho_prev_gen;
  
  // LO system
  std::vector<LA::MPI::SparseMatrix*> vec_lo_sys;
  std::vector<LA::MPI::Vector*> vec_lo_rhs;
  std::vector<LA::MPI::Vector*> vec_lo_fixed_rhs;
  std::vector<LA::MPI::Vector*> vec_lo_sflx;
  std::vector<LA::MPI::Vector*> vec_lo_sflx_old;
  std::vector<LA::MPI::Vector*> vec_lo_sflx_prev_gen;
  
  std::map<std::pair<unsigned int, unsigned int>, unsigned int> component_index;
  std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int> > inverse_component_index;
  std::map<std::pair<unsigned int, unsigned int>, unsigned int> reflective_direction_index;
  std::map<std::vector<unsigned int>, unsigned int> relative_position_to_id;

  std::string discretization;
  double total_angle;
  unsigned int n_azi;
  unsigned int n_total_ho_vars;
  unsigned int n_group;
  unsigned int n_material;
  int p_order;
  int global_refinements;
  double c_penalty;
  unsigned int n_dir;
  std::set<unsigned int> fissile_ids;
  std::vector<unsigned int> ncell_per_dir;
  std::vector<double> cell_size_all_dir;
  std::vector<double> axis_max_values;

  bool is_eigen_problem;
  bool do_nda;
  bool have_reflective_bc;
  bool is_explicit_reflective;
  std::unordered_map<int, bool> is_reflective_bc;
  std::unordered_map<int, bool> is_material_fissile;
  Table<1, unsigned int> material_id_table_1D;
  Table<2, unsigned int> material_id_table_2D;
  Table<3, unsigned int> material_id_table_3D;

  std::vector<Tensor<1, dim> > omega_i;
  std::vector<double> wi;
  std::vector<double> tensor_norms;
  
  double fission_source;
  double fission_source_prev_gen;

  const double err_k_tol;
  const double err_phi_tol;

  std::vector<std::vector<double> > all_sigt;
  std::vector<std::vector<double> > all_inv_sigt;
  std::vector<std::vector<double> > all_ksi;
  std::vector<std::vector<double> > all_nusigf;
  std::vector<std::vector<double> > all_q;
  std::vector<std::vector<double> > all_q_per_ster;

  std::vector<std::vector<std::vector<double> > > all_sigs;
  std::vector<std::vector<std::vector<double> > > all_sigs_per_ster;
  std::vector<std::vector<std::vector<double> > > all_ksi_nusigf;
  std::vector<std::vector<std::vector<double> > > all_ksi_nusigf_per_ster;
  std::vector<std::vector<std::vector<double> > > ho_scaled_fiss_transfer_per_ster;
  std::vector<std::vector<std::vector<double> > > lo_scaled_fiss_transfer;
 
  ConditionalOStream                        pcout;
  //TimerOutput                               computing_timer;
  
  std::vector<std_cxx11::shared_ptr<LA::MPI::PreconditionAMG> > pre_ho_amg;
};

template <int dim>
EP_SN<dim>::EP_SN (const ParameterHandler &prm)
:
mpi_communicator (MPI_COMM_WORLD),
triangulation (mpi_communicator,
               typename Triangulation<dim>::MeshSmoothing
               (Triangulation<dim>::smoothing_on_refinement |
                Triangulation<dim>::smoothing_on_coarsening)),
dof_handler (triangulation),
pi(dealii::numbers::PI),
discretization(prm.get("spatial discretization")),
n_azi(prm.get_integer("angular quadrature order")),
n_group(prm.get_integer("number of groups")),
n_material(prm.get_integer("number of materials")),
is_eigen_problem(prm.get_bool("do eigenvalue calculations")),
do_nda(prm.get_bool("do NDA")),
have_reflective_bc(prm.get_bool("have reflective BC")),
p_order(prm.get_integer("finite element polynomial degree")),
global_refinements(prm.get_integer("uniform refinements")),
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
}

template <int dim>
EP_SN<dim>::~EP_SN()
{
  dof_handler.clear();
  delete fe;
}

template <int dim>
void EP_SN<dim>::declare_parameters (ParameterHandler &prm0)
{
  // our final strategy is to declare all possible entries
  // and then ignore some of them suggested by Wolfgang Bangerth
  // from Colorado State on 05-10-2017
  // The following are the basic parameters we need to define a problem
  {
    prm0.declare_entry ("spatial discretization", "", Patterns::Selection("DFEM|DG|dfem|dg|CFEM|cfem|cg"), "USE DG or CG for spatial discretization");
    prm0.declare_entry ("angular quadrature order", "4", Patterns::Integer (), "Gauss-Chebyshev level-symmetric-like quadrature");
    prm0.declare_entry ("number of groups", "1", Patterns::Integer (), "Number of groups in MG calculations");
    prm0.declare_entry ("do eigenvalue calculations", "false", Patterns::Bool(), "Boolean to determine problem type");
    prm0.declare_entry ("do NDA", "false", Patterns::Bool(), "Boolean to determine NDA or not");
    prm0.declare_entry ("have reflective BC", "false", Patterns::Bool(), "");
    prm0.declare_entry ("reflective boundary names", "", Patterns::List (Patterns::Anything ()), "must be lower cases of xmin,xmax,ymin,ymax,zmin,zmax");
    prm0.declare_entry ("finite element polynomial degree", "1", Patterns::Integer(), "polynomial degree p for finite element");
    prm0.declare_entry ("uniform refinements", "0", Patterns::Integer(), "number of uniform refinements desired");
    prm0.declare_entry ("x, y, z max values of boundary locations", "", Patterns::List (Patterns::Double ()), "xmax, ymax, zmax of the boundaries, mins are zero");
    prm0.declare_entry ("number of cells for x, y, z directions", "", Patterns::List (Patterns::Integer ()), "Geotry is hyper rectangle defined by how many cells exist per direction");
    prm0.declare_entry("number of materials", "", Patterns::Integer (), "must be a positive integer");
    prm0.declare_entry("use explicit reflective boundary condition or not", "1", Patterns::Integer(), "");
    // prm0.declare_entry("material ID map", "", Patterns::List (Patterns::Integer ()), "Give material IDs for all blocks");
  }
  // FixIt: for current deal.II code, we don't consider reading mesh
  
  // Explanation: we brute-forcely declare as many entries as possible without read-in problem-definition
  // parameters. nmat and ngrp should both be large enough s.t. when reading starts, the real setting will
  // have entry-declaration
  unsigned int z_levels = 10;
  unsigned int y_levels = 100;
  prm0.enter_subsection ("Material ID map");
  {
    for (unsigned int z=0; z<z_levels; ++z)
      for (unsigned int y=0; y<y_levels; ++y)
      {
        std::ostringstream os;
        os << "z " << z << ", y " << y;
        prm0.declare_entry (os.str(), "", Patterns::List(Patterns::Integer()), "material IDs for a specific z level, y row");
      }
  }
  prm0.leave_subsection ();

  unsigned int nmat = 100;
  unsigned int ngrp = 10;
  prm0.enter_subsection ("sigma_t, group=1 to G");
  {
    for (unsigned int m=0; m<nmat; ++m)
    {
      std::ostringstream os;
      os << "material " << m + 1;
      prm0.declare_entry (os.str (), "", Patterns::List (Patterns::Double ()), "");
    }
  }
  prm0.leave_subsection ();
  
  for (unsigned int m=0; m<nmat; ++m)
  {
    std::ostringstream os;
    os << "sigma_s, material " << m + 1;
    prm0.enter_subsection (os.str());
    {
      for (unsigned int gin=0; gin<ngrp; ++gin)
      {
        std::ostringstream osm;
        osm << "g_in=" << gin + 1;
        prm0.declare_entry (osm.str(), "", Patterns::List(Patterns::Double()), "multigroup sigma_s");
      }
    }
    prm0.leave_subsection ();
  }

  prm0.enter_subsection ("one-group sigma_t");
  {
    for (unsigned int m=0; m<nmat; ++m)
    {
      std::ostringstream os;
      os << "material " << m + 1;
      prm0.declare_entry (os.str (), "1.0", Patterns::Double(), "");
    }
  }
  prm0.leave_subsection ();

  prm0.enter_subsection ("one-group sigma_s");
  {
    for (unsigned int m=0; m<nmat; ++m)
    {
      std::ostringstream os;
      os << "material " << m + 1;
      prm0.declare_entry (os.str (), "0.0", Patterns::Double(), "");
    }
  }
  prm0.leave_subsection ();

  prm0.enter_subsection ("one-group Q");
  {
    for (unsigned int m=0; m<nmat; ++m)
    {
      std::ostringstream os;
      os << "material " << m + 1;
      prm0.declare_entry (os.str (), "1.0", Patterns::Double(), "");
    }
  }
  prm0.leave_subsection ();
  
  prm0.enter_subsection ("Q, group=1 to G");
  {
    for (unsigned int m=0; m<nmat; ++m)
    {
      std::ostringstream os;
      os << "material " << m + 1;
      prm0.declare_entry (os.str (), "", Patterns::List (Patterns::Double ()), "");
    }
  }
  prm0.leave_subsection ();

  // the following is for eigen problems
  prm0.enter_subsection ("Fissile material IDs");
  {
    prm0.declare_entry ("fissile material IDs", "", Patterns::List (Patterns::Integer ()), "");
  }
  prm0.leave_subsection ();
  
  prm0.enter_subsection ("ksi, group=1 to G");
  {
    for (unsigned int m=0; m<nmat; ++m)
    {
      std::ostringstream os;
      os << "material " << m + 1;
      prm0.declare_entry(os.str(), "", Patterns::List(Patterns::Double()), "");
    } 
  }
  prm0.leave_subsection ();

  prm0.enter_subsection ("nu_sigf, group=1 to G");
  {
    for (unsigned int m=0; m<nmat; ++m)
    {
      std::ostringstream os;
      os << "material " << m + 1;
      prm0.declare_entry(os.str(), "", Patterns::List(Patterns::Double()), "");
    } 
  }
  prm0.leave_subsection ();
}

template <int dim>
void EP_SN<dim>::process_input (ParameterHandler &prm)
{
  process_problem_definition (prm);

  process_material_properties (prm);
}

template <int dim>
void EP_SN<dim>::process_problem_definition (ParameterHandler &prm)
{
  preprocess_reflective_bc (prm);

  process_coordinate_information (prm);
}

template <int dim>
void EP_SN<dim>::preprocess_reflective_bc (ParameterHandler &prm)
{
  if (have_reflective_bc)
  {
    is_explicit_reflective = prm.get_integer ("use explicit reflective boundary condition or not");
    // sorry, we need c++11
    std::map<std::string, unsigned int> bd_names_to_id {{"xmin",0},
                                                        {"xmax",1},
                                                        {"ymin",2},
                                                        {"ymax",3},
                                                        {"zmin",4},
                                                        {"zmax",5}};
    std::vector<std::string> strings = Utilities::split_string_list (prm.get ("reflective boundary names"));
    AssertThrow (strings.size()>0,
                 ExcMessage("reflective boundary names have to be entered"));
    std::set<int> tmp;
    for (unsigned int i=0; i<strings.size (); ++i)
    {
      AssertThrow(bd_names_to_id.find(strings[i])!=bd_names_to_id.end(),
                  ExcMessage("Invalid reflective boundary name: use xmin, xmax, etc."));
      tmp.insert (bd_names_to_id[strings[i]]);
    }
    auto it = tmp.begin ();
    std::ostringstream os;
    os << "No valid reflective boundary name for " << dim << "D";
    AssertThrow(*it<2*dim,
                ExcMessage(os.str()));
    for (unsigned int i=0; i<2*dim; ++i)
    {
      if (tmp.count (i))
        is_reflective_bc[i] = true;
      else
        is_reflective_bc[i] = false;
    }
  }
}

template <int dim>
void EP_SN<dim>::process_coordinate_information (ParameterHandler &prm)
{
  // max values for all axis
  {
    std::vector<std::string> strings = Utilities::split_string_list (prm.get ("x, y, z max values of boundary locations"));
    AssertThrow (strings.size()==dim,
                 ExcMessage("Number of axis max values must be the same as dimension"));
    for (unsigned int i=0; i<dim; ++i)
      axis_max_values.push_back (std::atof (strings[i].c_str()));
  }

  // read in number of cells and get cell sizes along axes
  {
    std::vector<std::string> strings = Utilities::split_string_list (prm.get ("number of cells for x, y, z directions"));
    AssertThrow (strings.size()==dim,
                 ExcMessage ("Entries for numbers of cells should be equal to dimension"));
    std::vector<unsigned int> cells_per_dir;
    std::vector<std::vector<double> > spacings;
    for (unsigned int d=0; d<dim; ++d)
    {  
      ncell_per_dir.push_back (std::atoi (strings[d].c_str ()));
      cell_size_all_dir.push_back (axis_max_values[d]/ncell_per_dir[d]);
    }
  }
}

template <int dim>
void EP_SN<dim>::initialize_relative_position_to_id_map (ParameterHandler &prm)
{
  prm.enter_subsection ("Material ID map");
  {
    for (unsigned int z=0; z<dim==3?ncell_per_dir[2]:1; ++z)
      for (unsigned int y=0; y<dim>=2?ncell_per_dir[1]:1; ++y)
      {
        std::ostringstream os;
        os << "z " << z << ", y " << y;
        std::vector<std::string> strings = Utilities::split_string_list (prm.get(os.str()));
        for (unsigned int x=0; x<ncell_per_dir[0]; ++x)
        {
          std::vector<unsigned int> tmp {x, y, z};
          relative_position_to_id[tmp] = std::atoi (strings[x].c_str());
        }
      }
  }
  prm.leave_subsection ();
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
void EP_SN<dim>::process_material_properties (ParameterHandler &prm)
{
  if (n_group>1)
  {
    // This block takes in sigts
    prm.enter_subsection ("sigma_t, group=1 to G");
    {
      for (unsigned int m=0; m<n_material; ++m)
      {
        std::ostringstream os;
        os << "material " << m + 1;
        std::vector<std::string> strings = Utilities::split_string_list (prm.get (os.str ()));
        AssertThrow (strings.size () == n_group,
                     ExcMessage ("n_group is not equal to group number of sigma_t"));
        std::vector<double> tmp, inv_tmp;
        for (unsigned int g=0; g<n_group; ++g)
        {
          tmp.push_back (std::atof (strings[g].c_str ()));
          inv_tmp.push_back (1.0/tmp[g]);
        }
        all_sigt.push_back (tmp);
        all_inv_sigt.push_back (inv_tmp);
      }
    }
    prm.leave_subsection ();
    
    // This block takes in scattering transfer matrices
    all_sigs.resize (n_material);
    all_sigs_per_ster.resize (n_material);
    for (unsigned int m=0; m<n_material; ++m)
    {
      std::ostringstream osm;
      osm << "sigma_s, material " << m + 1;
      std::vector<std::vector<double> >  tmp_sigs (n_group, std::vector<double>(n_group));
      std::vector<std::vector<double> >  tmp_sigs_per_ster (n_group, std::vector<double>(n_group));
      prm.enter_subsection (osm.str());
      {
        for (unsigned int gin=0; gin<n_group; ++gin)
        {
          std::ostringstream os;
          os << "g_in=" << gin + 1;
          std::vector<std::string> strings = Utilities::split_string_list (prm.get(os.str()));
          AssertThrow (strings.size()==n_group,
                       ExcMessage("sigma_s should have n_group entries per in group"));
          for (unsigned int g=0; g<n_group; ++g)
          {
            tmp_sigs[gin][g] = std::atof (strings[g].c_str());
            tmp_sigs_per_ster[gin][g] = std::atof (strings[g].c_str()) / (4.0 * pi);
          }
        }
      }
      prm.leave_subsection ();
      all_sigs[m] = tmp_sigs;
      all_sigs_per_ster[m] = tmp_sigs_per_ster;
    }
    
    if (!is_eigen_problem)
    {
      prm.enter_subsection ("Q, group=1 to G");
      {
        for (unsigned int m=0; m<n_material; ++m)
        {
          std::ostringstream os;
          os << "material " << m + 1;
          std::vector<std::string> strings = Utilities::split_string_list (prm.get (os.str ()));
          AssertThrow (strings.size () == n_group,
                       ExcMessage ("n_group is not equal to group number of Q"));
          std::vector<double> tmp_q, tmp_q_per_ster;
          for (unsigned int g=0; g<n_group; ++g)
          {
            tmp_q.push_back (std::atof (strings[g].c_str ()));
            tmp_q_per_ster.push_back (tmp_q[g] / (4.0 * pi));
          }
          all_q.push_back (tmp_q);
          all_q_per_ster.push_back (tmp_q_per_ster);
        }
      }
      prm.leave_subsection ();
    }
  }
  else
  {
    prm.enter_subsection ("one-group sigma_t");
    {
      std::vector<std::string> strings = Utilities::split_string_list (prm.get("values"));
      AssertThrow (strings.size()==n_material,
                   ExcMessage("One-group sigma_t should have N_material entries"));
      std::vector<double> tmp_sigt (n_material);
      for (unsigned int m=0; m<n_material; ++m)
      {
        // sorry, c++11 again.
        std::vector<double> tmp = {std::atof (strings[m].c_str())};
        all_sigt.push_back (tmp);
      }
    }
    prm.leave_subsection ();
    
    prm.enter_subsection ("one-group sigma_s");
    {
      std::vector<std::string> strings = Utilities::split_string_list (prm.get("values"));
      AssertThrow (strings.size()==n_material,
                   ExcMessage("One-group sigma_s should have N_material entries"));
      std::vector<double> tmp_sigt (n_material);
      for (unsigned int m=0; m<n_material; ++m)
      {
        // sorry, c++11 again.
        std::vector<std::vector<double> > tmp {std::vector<double> {std::atof(strings[m].c_str())}};
        all_sigs.push_back (tmp);
        std::vector<std::vector<double> > tmp_per_ster {std::vector<double> {std::atof(strings[m].c_str()) / (4.0*pi)}};
        all_sigs_per_ster.push_back (tmp_per_ster);
      }
    }
    prm.leave_subsection ();
    
    prm.enter_subsection ("one-group Q");
    {
      std::vector<std::string> strings = Utilities::split_string_list (prm.get("values"));
      AssertThrow (strings.size()==n_material,
                   ExcMessage("One-group sigma_s should have N_material entries"));
      std::vector<double> tmp_sigt (n_material);
      for (unsigned int m=0; m<n_material; ++m)
      {
        std::vector<double> tmp = {std::atof (strings[m].c_str())};
        all_q.push_back (tmp);
        std::vector<double> tmp_per_ster = {std::atof (strings[m].c_str()) / (4.0 * pi)};
        all_q_per_ster.push_back (tmp_per_ster);
      }
    }
    prm.leave_subsection ();
  }

  // This block is for eigenvalue problems
  if (is_eigen_problem)
  {
    process_eigen_material_properties (prm);
  } 
}

template <int dim>
void EP_SN<dim>::process_eigen_material_properties (ParameterHandler &prm)
{
  prm.enter_subsection ("Fissile material IDs");
  {
    std::ostringstream os;
    os << "fissile material IDs";
    std::vector<std::string> strings = Utilities::split_string_list (prm.get (os.str ()));
    AssertThrow (strings.size () > 0,
                 ExcMessage ("Fissile material IDs must be inserted for eigen problems"));
    // std::set<int> fissile_ids;
    for (unsigned int i=0; i<strings.size(); ++i)
      fissile_ids.insert (std::atoi (strings[i].c_str ()));
  }
  prm.leave_subsection ();

  for (unsigned int m=0; m<n_material; ++m)
  {
    if (fissile_ids.count(m))
      is_material_fissile[m] = true;
    else
      is_material_fissile[m] = false;
  }
  AssertThrow (!is_material_fissile.empty (),
               ExcMessage ("Please specify at least one valid ID for fissile materials"));

  if (n_group>1)
  {
    prm.enter_subsection ("ksi, group=1 to G");
    {
      for (unsigned int m=0; m<n_material;++m)
      {
        std::vector<double> tmp(n_group,0.0);
        if (is_material_fissile[m])
        {
          std::ostringstream os;
          os << "material " << m + 1;
          std::vector<std::string> strings = Utilities::split_string_list (prm.get (os.str ()));
          AssertThrow (strings.size () == n_group,
                       ExcMessage ("n_group is not equal to group number of ksi"));
          std::vector<double> tmp;
          for (unsigned int g=0; g<n_group; ++g)
            tmp[g] = std::atof (strings[g].c_str ());
        }
        all_ksi.push_back (tmp);
      }
    }
    prm.leave_subsection ();
    
    prm.enter_subsection ("nu_sigf, group=1 to G");
    {
      for (unsigned int m=0; m<n_material;++m)
      {
        std::vector<double> tmp(n_group,0.0);
        if (is_material_fissile[m])
        {
          std::ostringstream os;
          os << "material " << m + 1;
          std::vector<std::string> strings = Utilities::split_string_list (prm.get (os.str ()));
          AssertThrow (strings.size () == n_group,
                       ExcMessage ("n_group is not equal to group number of nusigf"));
          std::vector<double> tmp;
          for (unsigned int g=0; g<n_group; ++g)
            tmp[g] = std::atof (strings[g].c_str ());
        }
        all_nusigf.push_back (tmp);
      }
    }
    prm.leave_subsection ();
  }
  else
  {
    prm.enter_subsection ("one-group ksi");
    {
      std::vector<std::string> strings = Utilities::split_string_list (prm.get("values"));
      AssertThrow (strings.size()==n_material,
                   ExcMessage("One-group ksi should have N_material entries"));
      std::vector<double> tmp_sigt (n_material);
      for (unsigned int m=0; m<n_material; ++m)
      {
        // sorry, c++11 again.
        std::vector<double> tmp {is_material_fissile[m] ? std::atof (strings[m].c_str()) : 0.0};
        all_sigt.push_back (tmp);
      }
    }
    prm.leave_subsection ();
    
    prm.enter_subsection ("one-group nu_sigf");
    {
      std::vector<std::string> strings = Utilities::split_string_list (prm.get("values"));
      AssertThrow (strings.size()==n_material,
                   ExcMessage("One-group sigma_t should have N_material entries"));
      std::vector<double> tmp_sigt (n_material);
      for (unsigned int m=0; m<n_material; ++m)
      {
        std::vector<double> tmp {is_material_fissile[m] ? std::atof (strings[m].c_str()) : 0.0};
        all_sigt.push_back (tmp);
      }
    }
    prm.leave_subsection ();
  }

  for (unsigned int m=0; m<n_material; ++m)
  {
    std::vector<std::vector<double> >  tmp (n_group, std::vector<double>(n_group));
    std::vector<std::vector<double> >  tmp_per_ster (n_group, std::vector<double>(n_group));
    if (is_material_fissile[m])
      for (unsigned int gin=0; gin<n_group; ++gin)
        for (unsigned int g=0; g<n_group; ++g)
        {
          tmp[gin][g] = all_ksi[m][g] * all_nusigf[m][gin];
          tmp_per_ster[gin][g] = tmp[gin][g] / (4.0 * pi);
        }
    all_ksi_nusigf.push_back (tmp_per_ster);
    all_ksi_nusigf_per_ster.push_back (tmp_per_ster);
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
void EP_SN<dim>::initialize_component_index ()
{
  // initialize the map from (direction, group) to component indices
  unsigned int ind = 0;
  for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
    for (unsigned int g=0; g<n_group; ++g)
    {
      std::pair<unsigned int, unsigned int> key (i_dir, g);
      component_index.insert (std::make_pair (key, ind));
      inverse_component_index[ind] = key;
      ind += 1;
    }
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
void EP_SN<dim>::initialize_ref_bc_index ()
{
  AssertThrow (dim>1,
               ExcMessage("1D cases are not implemented for now."));
  std::vector<Tensor<1, dim> > boundary_normal_vectors;
  boundary_normal_vectors.size (2*dim);
  if (dim==2)
  {
    boundary_normal_vectors[0][0] = -1.0;
    boundary_normal_vectors[1][0] = 1.0;
    boundary_normal_vectors[2][1] = -1.0;
    boundary_normal_vectors[3][1] = 1.0;
  }
  if (dim==3)
  {
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
    
    vec_ho_rhs.push_back (new LA::MPI::Vector);
    vec_ho_sflx.push_back (new LA::MPI::Vector);
    vec_ho_sflx_old.push_back (new LA::MPI::Vector);
    vec_ho_fixed_rhs.push_back (new LA::MPI::Vector);

    for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
    {
      vec_ho_sys.push_back (new LA::MPI::SparseMatrix);
      vec_aflx.push_back (new LA::MPI::Vector);
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
            cell->face(fn)->set_boundary_id (0);
          
          // right boundary
          if (std::fabs(ct[0]-axis_max_values[0])<1.0e-14)
            cell->face(fn)->set_boundary_id (1);
          
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
  assemble_ho_volume_boundary ();
  
  if (discretization=="DFEM")
    assemble_ho_interface ();
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
  
  // volumetric matrix
  std::vector<FullMatrix<double> > cell_matrix(n_total_ho_vars, FullMatrix<double>(dofs_per_cell, dofs_per_cell));
  FullMatrix<double> cell_matrix_collision(dofs_per_cell, dofs_per_cell);
  
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  
  // volumetric pre-assembly matrices
  std::vector<FullMatrix<double> > mass_at_qp (n_q, FullMatrix<double>(dofs_per_cell, dofs_per_cell));
  std::vector<std::vector<FullMatrix<double> > > stiffness_at_qp (n_q, 
      std::vector<FullMatrix<double> > (n_dir, FullMatrix<double> (dofs_per_cell, dofs_per_cell)));
  
  bool pre_assemble_cell_finished = false;
  
  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (typename DoFHandler<dim>::active_cell_iterator cell=dof_handler.begin_active(); 
      cell!=dof_handler.end(); ++cell)
  {
    if (cell->is_locally_owned())
    {
      fv.reinit(cell);
      cell->get_dof_indices (local_dof_indices);
      std::vector<double> local_sigt = all_sigt[cell->material_id()];
      std::vector<double> local_inv_sigt = all_inv_sigt[cell->material_id()];
      std::vector<FullMatrix<double> > local_matrices (n_total_ho_vars, FullMatrix<double> (dofs_per_cell, dofs_per_cell));
      // FixIt: a more proper definition for h_cell
      
      if (!pre_assemble_cell_finished)
      {
        for (unsigned int qi=0; qi<n_q; ++qi)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              mass_at_qp[qi](i,j) += (fv.shape_value(i,qi) *
                                      fv.shape_value(j,qi));
        
        
        for (unsigned int qi=0; qi<n_q; ++qi)
          for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              for (unsigned int j=0; j<dofs_per_cell; ++j)
                stiffness_at_qp[qi][i_dir](i,j) += ((fv.shape_grad(i,qi) *
                                                     omega_i[i_dir]) 
                                                    *
                                                    (fv.shape_grad(j,qi) *
                                                     omega_i[i_dir]));
        
        pre_assemble_cell_finished = true;
      }
      
      // Use pre-assembled matrix components in reference cell to do assembly in real matrix
      FullMatrix<double> unscaled_mass(dofs_per_cell, dofs_per_cell);
      std::vector<FullMatrix<double> > unscaled_stiffness(n_dir, FullMatrix<double>(dofs_per_cell, dofs_per_cell));
      
      for (unsigned int qi=0; qi<n_q; ++qi)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++i)
            unscaled_mass(i,j) += mass_at_qp[qi](i,j) * fv.JxW(qi);
        
      for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
        for (unsigned int qi=0; qi<n_q; ++qi)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              unscaled_stiffness[i_dir](i,j) += stiffness_at_qp[qi][i_dir](i,j) * fv.JxW(qi);
          
      for (unsigned int g=0; g<n_group; ++g)
      {
        FullMatrix<double> scaled_mass (dofs_per_cell, dofs_per_cell);
        for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
        {
          if (i_dir==0)
            scaled_mass.add (local_sigt[g], unscaled_mass);
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
            for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
            {
              double absndo;
              // Note: we assume the face is not curvilinear
              absndo = std::fabs(vec_n * omega_i[i_dir]);
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
          }
          else
          {
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
                        local_matrices[ind](i,j) += -(absndo *
                                                      fvf.shape_value(i,qi) *
                                                      (omega_i[r_dir] * fvf.shape_grad(j,qi)) /
                                                      local_sigt[g] *
                                                      fvf.JxW(qi));
                }// g
              }// i_dir
          }// is_reflective_bc
        }// boundary faces for robin boundaries
      }// face
      
      for (unsigned int k=0; k<n_total_ho_vars; ++k)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            vec_ho_sys[k]->add (local_dof_indices[i],
                                local_dof_indices[j],
                                local_matrices[k](i, j));
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
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              for (unsigned int j=0; j<dofs_per_cell; ++j)
              {
                vec_ho_sys[k]->add (local_dof_indices[i],
                                    local_dof_indices[j],
                                    all_real_vp_up[k](i,j));
                
                vec_ho_sys[k]->add (local_dof_indices[i],
                                    neigh_dof_indices[j],
                                    all_real_vp_un[k](i,j));
                
                vec_ho_sys[k]->add (neigh_dof_indices[i],
                                    local_dof_indices[j],
                                    all_real_vn_up[k](i,j));
                
                vec_ho_sys[k]->add (neigh_dof_indices[i],
                                    neigh_dof_indices[j],
                                    all_real_vn_un[k](i,j));
              }// j
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
      tmp[0] = mu_quad.point(i)[0] * 2.0 - 1.0;
      omega_i.push_back (tmp);
      wi.push_back (mu_quad.weight(i));
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
      unsigned int level_angle_num = 4 * (n_azi - i);
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
      unsigned int level_angle_num = i < n_azi / 2 ? 4 * (i + 1) : 4 * (n_azi - i);
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
  n_total_ho_vars = n_dir * n_group;
  quadr.close();
  
  /*
  for (unsigned int i=0; i<n_total_ho_vars; ++i)
  {
    FEValuesExtractors::Scalar tmp(i);
    comp.push_back(tmp);
  }
   */
  
  // estimate tensor norm to do penalty method
  for (unsigned int i=0; i<n_dir; ++i)
  {
    Tensor<2, dim> tensor_tmp = outer_product(omega_i[i], omega_i[i]);
    tensor_norms.push_back(tensor_tmp.norm());
  }
}

template <int dim>
void EP_SN<dim>::initialize_ho_preconditioners ()
{
  //TimerOutput::Scope t (computing_timer, "HO preconditioner initialization");
  pre_ho_amg.resize (n_total_ho_vars);
  for (unsigned int i=0; i<n_total_ho_vars; ++i)
  {
    pre_ho_amg[i].reset ();
    pre_ho_amg[i] = (std_cxx11::shared_ptr<LA::MPI::PreconditionAMG> (new LA::MPI::PreconditionAMG));
    LA::MPI::PreconditionAMG::AdditionalData data;
    if (!is_explicit_reflective)
      data.symmetric_operator = true;
    else
      data.symmetric_operator = false;
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
                                  vec_ho_rhs[i]->l2_norm() * 1e-12);

#ifdef USE_PETSC_LA
    if (!is_explicit_reflective)
    {
      LA::SolverCG solver (solver_control, mpi_communicator);
      solver.solve (*(vec_ho_sys)[i],
                    *(vec_aflx)[i],
                    *(vec_ho_rhs)[i],
                    *(pre_ho_amg)[i]);
    }
    else
    {
      LA::SolverBicgstab solver (solver_control, mpi_communicator);
      solver.solve (*(vec_ho_sys)[i],
                    *(vec_aflx)[i],
                    *(vec_ho_rhs)[i],
                    *(pre_ho_amg)[i]);
    }
#else
    // LA::SolverCG solver(solver_control);
#endif
  }  
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
    for (unsigned int g=0; g<n_group; ++g)
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
    for (unsigned int g=0; g<n_group; ++g)
    {
      *(vec_ho_sflx_old)[g] = *(vec_ho_sflx)[g];
      *(vec_ho_sflx)[g] = 0;
      for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
        vec_ho_sflx[g]->add(wi[i_dir], *(vec_aflx)[get_component_index(i_dir, g)]);
    }
}

template <int dim>
void EP_SN<dim>::generate_ho_source ()
{
  const QGauss<dim>  q_rule(p_order+1);
  const QGauss<dim-1>  qf_rule(p_order+1);
  
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
      std::vector<Vector<double> > vec_cell_rhs(n_group, Vector<double> (dofs_per_cell));
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
              if (all_sigs_per_ster[material_id][gin][g]>1.0e-13)
                vec_cell_rhs[g](i) += test_func_jxw * all_sigs_per_ster[material_id][gin][g];
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
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          vec_ho_rhs[k]->add(local_dof_indices,
                             vec_cell_rhs[k]);
    }// local cells
  }

  for (unsigned int i=0; i<n_total_ho_vars; ++i)
    vec_ho_rhs[i]->compress (VectorOperation::add);
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
  std::vector<Vector<double> > vec_cell_rhs(n_group, Vector<double> (dofs_per_cell));
 
  for (typename DoFHandler<dim>::active_cell_iterator cell = triangulation.begin_active(); cell!= triangulation.end(); ++cell)
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
            fv.get_function_values (vec_ho_sflx[g], local_ho_sflxes[g]);

          for (unsigned int qi=0; qi<n_q; ++qi)
          {
            double test_func_jxw;
            for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              test_func_jxw = fv.shape_value (i, qi) * fv.JxW (qi);
              for (unsigned int g=0; g<n_group; ++g)
                for (unsigned int gin=0; g<n_group; ++gin)
                  vec_cell_rhs[g](i) += test_func_jxw * ho_scaled_fiss_transfer_per_ster[material_id][gin][g];
            }
          }

          for (unsigned int g=0; g<n_group; ++g)
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              vec_ho_fixed_rhs[g](local_dof_indices[i]) += vec_cell_rhs[g](i);
        }
      }
      else
      {
        if (*(std::max_element(all_q_per_ster[material_id]))>1.0e-13)
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

  initialize_ho_preconditioners ();

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
  //TimerOutput::Scope t (computing_timer, sec_name);
  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);

  for (unsigned int g=0; g<n_group; ++g)
  {
    std::ostringstream os;
    os << "ho-phi g=" << g;
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
void EP_SN<dim>::run (ParameterHandler &prm)
{
  process_input (prm);
  
  generate_globally_refined_grid ();
  
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

int main(int argc, char *argv[])
{
  try
  {
    using namespace dealii;
 
    int dimension;
    if (argc!=3)
    {
      std::cerr << "Call the program as ./dg-ep-proto input_file_name dimension" << std::endl;
      return 1;
    }
    else
    {
      std::stringstream convert(argv[2]);
      convert >> dimension;
      assert (dimension == 2 || dimension == 3);
    }
    ParameterHandler prm;
    switch (dimension)
    {
      case 2:
      {
        EP_SN<2>::declare_parameters(prm);
        prm.read_input(argv[1]);
        Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
        EP_SN<2> dg_ep (prm);
        dg_ep.run (prm);
        break;
      }
        
      case 3:
      {
        EP_SN<3>::declare_parameters(prm);
        prm.read_input(argv[1]);
        Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
        EP_SN<3> dg_ep (prm);
        dg_ep.run (prm);
        break;
      }
        
      default:
        AssertThrow (dimension>1,
                     ExcMessage("1D is not implemented for now."));
    }
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
