#ifndef __preconditioner_solver_h__
#define __preconditioner_solver_h__

#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/base/parameter_handler.h>

#include <vector>
#include <string>

using namespace dealii;

//! The class provides linear solving functionalities.
/*!
 This class is a wrapper class aiming to simplifies linear solving process. deal.II
 algebraic preconditioner wrappers and solver wrappers will be called to utilize
 PETSc linear solving functionalities in parallel. This class only provide two
 functions interfaces outside: initialize_preconditioners and linear_algebraic_solve.
 
 About how to use this class:
 
 (1) Everytime a new equation is assembled, call initialize_preconditioners to
     either initilize preconditioners if iterative solvers are used or initilize
     factorization for direct solver based on MUMPS. Note that initialization only
     needs to be done once per equation.
 
 (2) Call linear_algebra_solve whenever needed.
 
 \author Weixiong Zheng
 \date: 2017/10
 */
class PreconditionerSolver
{
public:
  /*!
   Class constructor.
   
   \param prm A ParameterHandler object containing all user defined parameters.
   \param equation_name A string describing the name of the target equation.
   \param n_total_vars A integer for the total number components for target equation.
   */
  PreconditionerSolver (const ParameterHandler &prm,
                        std::string equation_name,
                        unsigned int& n_total_vars);
  
  //! Class destructors.
  ~PreconditionerSolver ();
  
  /*!
   This function initilizes preconditioners for an equation. The main 
   functionalities include one of the following two points:
   
   (1) For iterative solvers, initializing preconditioners and store them in 
       shared_ptr's.
   
   (2) For MUMPS direct solver, initializing factorizations and store them in
       shared_ptr's.
   
   \param sys_mats A vector of pointers to system matrices of the target equation.
   \param sys_rhses A vector of pointers to system right-hand-side vectors of the
          target equation.
   */
  void initialize_preconditioners
  (std::vector<PETScWrappers::MPI::SparseMatrix*> &sys_mats,
   std::vector<PETScWrappers::MPI::Vector*> &sys_rhses);
  
  /*!
   This function provide functionality of performing linear algebraic solve for 
   a specific component of the target equation.
   
   \param sys_mats A vector of pointers to system matrices of the target equation.
   \param sys_flxes A vector of pointers to system solutions of the target equation
   \param sys_rhses A vector of pointers to system right-hand-side vectors of the
          target equation.
   \param i An integer specifying the component index of interest.
   \return Void.
   */
  void linear_algebra_solve
  (std::vector<PETScWrappers::MPI::SparseMatrix*> &sys_mats,
   std::vector<PETScWrappers::MPI::Vector*> &sys_flxes,
   std::vector<PETScWrappers::MPI::Vector*> &sys_rhses,
   unsigned int &i/*component number*/);

private:
  const unsigned int n_total_vars;//!< Total number of variables in the target system.
  const std::string equation_name;//!< Name of the target equation.
  
  bool have_reflective_bc;//!< Booolean to determine if there's any reflective boundary.
  double ssor_omega;//!< The relaxation factor if block SSOR is used as preconditioner.
  
  std::string linear_solver_name;//!< Linear solver name for current equation.
  std::string preconditioner_name;//!< Preconditioner name for current equation.
  
  //! A vector of boolean to determine if direct solver of a component is initilized.
  std::vector<bool> direct_init;
  
  //! A vector of integers showing number of iterations in linear solves.
  std::vector<unsigned int> linear_iters;
  
  //! pointer of SolverControl object.
  /*!
   Control iterative algebraic linear solvers to determine convergence. Can be used
   to show number of linear iterations in the algebraic solve, linear solver residual
   after each iteration etc.
   */
  std_cxx11::shared_ptr<SolverControl> cn;
  
  //! A vector of pointers of BoomerAMG preconditioner.
  std::vector<std_cxx11::shared_ptr<PETScWrappers::PreconditionBoomerAMG> > pre_amg;
  
  //! A vector of pointers of block Jacobi preconditioner.
  std::vector<std_cxx11::shared_ptr<PETScWrappers::PreconditionBlockJacobi> > pre_bjacobi;
  
  //! A vector of pointers of ParaSails preconditioner.
  std::vector<std_cxx11::shared_ptr<PETScWrappers::PreconditionParaSails> > pre_parasails;
  
  //! A vector of pointers of Jacobi preconditioner.
  std::vector<std_cxx11::shared_ptr<PETScWrappers::PreconditionJacobi> > pre_jacobi;
  
  //! A vector of pointers of Eisenstat (block SSOR) preconditioner.
  std::vector<std_cxx11::shared_ptr<PETScWrappers::PreconditionEisenstat> > pre_eisenstat;
  
  //! A vector of pointers of MUMPS solver objects.
  std::vector<std_cxx11::shared_ptr<PETScWrappers::SparseDirectMUMPS> > direct;
};

#endif //__preconditioner_solver_h__
