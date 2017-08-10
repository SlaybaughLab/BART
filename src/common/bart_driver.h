#ifndef __transport_base_h__
#define __transport_base_h__
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

#include "../iteration/Iterations.h"
#include "../common/problem_definition.h"
#include "../common/preconditioner_solver.h"
#include "../mesh/mesh_generator.h"
#include "../material/material_properties.h"
#include "../aqdata/aq_base.h"

using namespace dealii;

template <int dim>
class BartDriver
{
public:
  BartDriver (ParameterHandler &prm);// : ProblemDefinition<dim> (prm){}
  virtual ~BartDriver ();
  
  void run ();
private:
  void setup_system ();
  void report_system ();
  void initialize_dealii_objects ();
  void prepare_correction_aflx ();
  void output_results () const;
  void initialize_aq (ParameterHandler &prm);
  
  const ParameterHandler &prm;
  
  Iterations<dim> itr_cls;
  std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr;
  std_cxx11::shared_ptr<MaterialProperties> mat_ptr;
  std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr;
  std_cxx11::shared_ptr<PreconditionerSolver> sol_ptr;
  
  std::string transport_model_name;
  std::string ho_linear_solver_name;
  std::string ho_preconditioner_name;
  std::string discretization;
  std::string namebase;
  std::string aq_name;

  FE_Poly<TensorProductPolynomials<dim>,dim,dim>* fe;
  
  parallel::distributed::Triangulation<dim> triangulation;
  
  DoFHandler<dim> dof_handler;
  
  IndexSet local_dofs;
  IndexSet relevant_dofs;
  
  double keff;
  
  bool is_eigen_problem;
  bool do_nda;
  bool have_reflective_bc;
  
  
  unsigned int n_dir;
  unsigned int n_azi;
  unsigned int n_total_ho_vars;
  unsigned int n_group;
  unsigned int n_material;
  unsigned int p_order;
  unsigned int global_refinements;
  
  ConstraintMatrix constraints;
};

#endif	// define  __transport_base_h__
