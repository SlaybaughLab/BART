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

#include <iostream>
#include <utility>
#include <vector>

#include "eigen_base.h"
#include "mg_base.h"
#include "../common/problem_definition.h"
#include "../common/preconditioner_solver.h"
#include "../mesh/mesh_generator.h"
#include "../material/material_properties.h"
#include "../aqdata/aq_base.h"
#include "../equation/equation_base.h"

using namespace dealii;

template <int dim>
class Iterations
{
public:
  Iterations (const ParameterHandler &prm,
              const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
              const std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr,
              const std_cxx11::shared_ptr<MaterialProperties> mat_ptr);
  ~Iterations ();
  
  void solve_problems (std::vector<Vector<double> > &sflxes_proc);
  
  void initialize_cell_iterators_this_proc
  (const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
   const DoFHandler<dim> &dof_handler);
  
  void initialize_assembly_related_objects
  (FE_Poly<TensorProductPolynomials<dim>,dim,dim>* fe);
  
  void initialize_system_matrices_vectors
  (DynamicSparsityPattern &dsp,
   IndexSet &local_dofs,
   std::vector<Vector<double> > &sflxes_proc);
  
  void get_keff (double &keff);

private:
  const std::string transport_name;
  
  std_cxx11::shared_ptr<EigenBase<dim> > eig_ptr;
  std_cxx11::shared_ptr<MGBase<dim> > mg_ptr;
  
  std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > equ_ptrs;
  
  double keff;
  bool is_eigen_problem;
};

#endif	// define  __iterations_h__
