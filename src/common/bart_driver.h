#ifndef __bart_driver_h__
#define __bart_driver_h__
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

#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../mesh/mesh_generator.h"
#include "../material/material_properties.h"
#include "../aqdata/aq_base.h"
#include "../equation/equation_base.h"
#include "../iteration/mg_base.h"
#include "../iteration/ig_base.h"
#include "../iteration/eigen_base.h"
#include "../iteration/iterations.h"

using namespace dealii;

template <int dim>
class BartDriver
{
public:
  BartDriver (ParameterHandler &prm);// : ProblemDefinition<dim> (prm){}
  ~BartDriver ();
  
  void run ();
private:
  void build_basis (ParameterHandler &prm);
  void setup_system ();
  void report_system ();
  void initialize_dealii_objects ();
  void output_results () const;
  
  void build_equation
  (std_cxx11::shared_ptr<EquationBase<dim> > equ_ptr,
   std::string equation_name,
   const ParameterHandler &prm,
   const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
   const std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr,
   const std_cxx11::shared_ptr<MaterialProperties> mat_ptr);
  
  void build_aq_model
  (std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr, ParameterHandler &prm);
  
  /** \brief Function used to build pointer to instance of InGroupBase's derived class
   */
  void build_eigen_iterations
  (std_cxx11::shared_ptr<EigenBase<dim> > eig_ptr, const ParameterHandler &prm);
  
  /** \brief Function used to build pointer to instance of MGBase's derived class
   */
  void build_mg_iterations
  (std_cxx11::shared_ptr<MGBase<dim> > mg_ptr, const ParameterHandler &prm);
  
  /** \brief Function used to build pointer to instance of InGroupBase's derived class
   */
  void build_ig_iterations
  (std_cxx11::shared_ptr<IGBase<dim> > ig_ptr, const ParameterHandler &prm);
  
  void initialize_cell_iterators_this_proc
  (const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
   const DoFHandler<dim> &dof_handler);
  
  void initialize_assembly_related_objects
  (FE_Poly<TensorProductPolynomials<dim>,dim,dim>* fe);
  
  void initialize_system_matrices_vectors
  (DynamicSparsityPattern &dsp,
   IndexSet &local_dofs,
   std::vector<Vector<double> > &sflxes_proc);
  
  std_cxx11::shared_ptr<Iterations<dim> > itr_ptr;
  std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr;
  std_cxx11::shared_ptr<MaterialProperties> mat_ptr;
  std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr;
  std_cxx11::shared_ptr<EigenBase<dim> > eig_ptr;
  std_cxx11::shared_ptr<MGBase<dim> > mg_ptr;
  std_cxx11::shared_ptr<IGBase<dim> > ig_ptr;
  std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > equ_ptrs;
  
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
  
  std::vector<Vector<double> > sflxes_proc;
  
  ConstraintMatrix constraints;
  ConditionalOStream pcout;
};

#endif	// define  __bart_driver_h__
