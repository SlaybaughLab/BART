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
#include "ig_base.h"
#include "../common/problem_definition.h"
#include "../common/preconditioner_solver.h"
#include "../equation/equation_base.h"


using namespace dealii;

//! This class provide an interface between BartDriver and specific iteration schemes.
/*!
 This class implements an interface between BartDriver<dim> and iteration schemes.
 The main purpose is such that one would not need to develop different functions
 according to different problem types.
 */
template <int dim>
class Iterations
{
public:
  /*!
   Class constructor.
   
   \param prm const ParameterHandler object.
   */
  Iterations (const ParameterHandler &prm);
  
  //! Class destructor.
  ~Iterations ();
  
  /*!
   A generic function to solve the problem for both eigenvalue and fixed-source
   problems. Inside this function, solving methods will be specified according to
   Iterations<dim>::is_eigen_problem the boolean.
   
   \note In fixed-source problems, eig_ptr points to a null pointer and will do 
   nothing.
   
   \param equ_ptrs A vector of shared_ptr's of EquationBase objects.
   \param ig_ptr A shared_ptr of downcast IGBase object.
   \param mg_ptr A shared_ptr of downcast MGBase object.
   \param eig_ptr A shared_ptr of downcast EigenBase object.
   \return Void.
   */
  void solve_problems
  (std::vector<Vector<double> > &sflxes_proc,
   std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs,
   std_cxx11::shared_ptr<IGBase<dim> > ig_ptr,
   std_cxx11::shared_ptr<MGBase<dim> > mg_ptr,
   std_cxx11::shared_ptr<EigenBase<dim> > eig_ptr);
  
  /*!
   A function to retrieve value of keff.
   
   \todo Change this function to return fashion to keep consistency of get_xxx type
   functions.
   
   \param keff A double living in BartDriver<dim> to be modified in-place.
   \return Void.
   */
  void get_keff (double &keff);

private:
  const std::string transport_name;//!< HO equation name.
  
  double keff;//!< keff living in current class.
  bool is_eigen_problem;//!< Boolean to determine if the problem is eigenvalue problem.
};

#endif	// define  __iterations_h__
