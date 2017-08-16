#ifndef __iterations_h__
#define __iterations_h__
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/numerics/vector_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_poly.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <unordered_map>
#include <utility>
#include <vector>

#include "../common/problem_definition.h"
#include "../common/preconditioner_solver.h"
#include "../mesh/mesh_generator.h"
#include "../material/material_properties.h"
#include "../aqdata/aq_base.h"

using namespace dealii;

template <int dim>
class Iterations
{
public:
  Iterations (ParameterHandler &prm);// : ProblemDefinition<dim> (prm){}
  virtual ~Iterations ();
  
  void do_iterations ();
  
  void initialize_system_matrices_vectors
  (SparsityPatternType &dsp, IndexSet &local_dofs);
  
private:
  void power_iteration ();
  void initialize_fiss_process ();
  void update_ho_moments_in_fiss ();
  void update_fiss_source_keff ();
  void source_iteration ();
  void scale_fiss_transfer_matrices ();
  void renormalize_sflx (std::vector<LA::MPI::Vector*> &target_sflxes);
  void NDA_PI ();
  void NDA_SI ();
  void get_flux_this_proc (std::vector<Vector<double> > &sflxes_proc);
  void get_keff (double &keff);
  
  double estimate_phi_diff (std::vector<LA::MPI::Vector*> &phis_newer,
                            std::vector<LA::MPI::Vector*> &phis_older);
  
  std_cxx11::shared_ptr<MaterialProperties> mat_ptr;
  std_cxx11::shared_ptr<PreconditionerSolver> sol_ptr;
  
  std::string namebase;
  std::string aq_name;
  
protected:
  std::vector<typename DoFHandler<dim>::active_cell_iterator> local_cells;
  std::vector<typename DoFHandler<dim>::active_cell_iterator> ref_bd_cells;
  std::vector<bool> is_cell_at_bd;
  std::vector<bool> is_cell_at_ref_bd;

  const double err_k_tol;
  const double err_phi_tol;
  const double err_phi_eigen_tol;
  
  double keff;
  double keff_prev_gen;
  double fission_source;
  double fission_source_prev_gen;
  
  bool is_eigen_problem;
  bool do_nda;
  
  unsigned int n_dir;
  unsigned int n_azi;
  unsigned int n_total_ho_vars;
  unsigned int n_group;
  unsigned int n_material;
  unsigned int p_order;
  
  std::vector<unsigned int> linear_iters;
  
  // HO system
  \std::vector<LA::MPI::Vector*> vec_aflx;
  std::vector<LA::MPI::Vector*> vec_ho_rhs;
  std::vector<LA::MPI::Vector*> vec_ho_fixed_rhs;
  std::vector<LA::MPI::Vector*> vec_ho_sflx_old;
  std::vector<LA::MPI::Vector*> vec_ho_sflx_prev_gen;
  
  // LO system
  std::vector<LA::MPI::Vector*> vec_lo_rhs;
  std::vector<LA::MPI::Vector*> vec_lo_fixed_rhs;
  std::vector<LA::MPI::Vector*> vec_lo_sflx;
  std::vector<LA::MPI::Vector*> vec_lo_sflx_old;
  std::vector<LA::MPI::Vector*> vec_lo_sflx_prev_gen;
  
  std::vector<Vector<double> > sflx_proc;
  std::vector<Vector<double> > sflx_proc_prev_gen;
  std::vector<Vector<double> > lo_sflx_proc;
};

#endif	// define  __iterations_h__
