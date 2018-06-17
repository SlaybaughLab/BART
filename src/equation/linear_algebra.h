#ifndef BART_SRC_EQUATION_LINEAR_ALGEBRA_H_
#define BART_SRC_EQUATION_LINEAR_ALGEBRA_H_

//! The class provides linear solving functionalities.
/*!
 This class is a wrapper class aiming to simplifies linear solving process. deal.II
 algebraic preconditioner wrappers and solver wrappers will be called to utilize
 PETSc linear solving functionalities in parallel. This class only provide two
 functions interfaces outside: InitPrecond and LinAlgSolve.

 About how to use this class:

 (1) Everytime a new equation is assembled, call InitPrecond to either initilize
     preconditioners if iterative solvers are used or initilize factorization
     for direct solver based on MUMPS. Note that initialization only needs to be
     done once per equation.

 (2) Call LinAlgSolve whenever needed.

 For details of PETSc wrappers of linear algebra functions/objects, refer to
 <a href="https://www.dealii.org/8.5.0/doxygen/deal.II/group__dealii::PETScWrappers.html"
 style="color:blue"><b>dealii::PETScWrappers</b></a>.

 \author Weixiong Zheng
 \date 2017, 2018
 */
class LinearAlgebra {
 public:
  /*!
   Class constructor.

   \param prm A ParameterHandler object containing all user defined parameters.
   \param equation_name A string describing the name of the target equation.
   \param n_total_vars A integer for the total number components for target equation.
   */
  LinearAlgebra (const dealii::ParameterHandler &prm,
      const std::string &equ_name,
      const int &n_total_vars);

  //! Class destructors.
  ~LinearAlgebra ();

  /*!
   This function initilizes preconditioners for an equation. The main
   functionalities include one of the following two points:

   (1) For iterative solvers, initializing preconditioners and store them in
       unique_ptr's.

   (2) For MUMPS direct solver, initializing factorizations and store them in
       unique_ptr's.

   \param sys_mats Hash table containing pointers to system matrices of the target equation.
   \param sys_rhses Hash table containing pointers to system right-hand-side vectors of the
          target equation.
   */
  void InitPrecond (std::vector<dealii::PETScWrappers::MPI::SparseMatrix*> &sys_mats,
      std::vector<dealii::PETScWrappers::MPI::Vector*> &sys_rhses);

  /*!
   This function provide functionality of performing linear algebraic solve for
   a specific component of the target equation.

   \param sys_mats Hash table containing pointers to system matrices of the target equation.
   \param sys_flxes Hash table containing pointers to system solutions of the target equation
   \param sys_rhses Hash table containing pointers to system right-hand-side vectors of the
          target equation.
   \param i An integer specifying the component index of interest.
   \param constraints Constraint matrix for local refinement.
   \return Void.
   */
  void LinAlgSolve (
      std::vector<dealii::PETScWrappers::MPI::SparseMatrix*> &sys_mats,
      std::vector<dealii::PETScWrappers::MPI::Vector*> &sys_flxes,
      std::vector<dealii::PETScWrappers::MPI::Vector*> &sys_rhses,
      std::vector<dealii::ConstraintMatrix*> &constraints,
      const int &i);

private:
  const int n_total_vars_;//!< Total number of variables in the target system.

  const bool have_reflective_bc_;//!< Boolean to determine if there's any reflective boundary

  const std::string equ_name_;//!< Name of the target equation.

  double ssor_omega_;//!< The relaxation factor if block SSOR is used as preconditioner.

  std::string linear_solver_name_;//!< Linear solver name for current equation.
  std::string preconditioner_name_;//!< Preconditioner name for current equation.

  /*!
   A Hash table having boolean to determine if direct solver of a component
   is initilized.
   */
  std::unordered_map<int, bool> direct_init_;

  //! Hash table containing integers showing number of iterations in linear solves.
  std::unordered_map<int, int> linear_iters_;

  //! pointer of <a href="https://www.dealii.org/8.4.1/doxygen/deal.II/classSolverControl.html" style="color:blue"><b>SolverControl</b></a> object.
  /*!
   Control iterative algebraic linear solvers to determine convergence. Can be used
   to show number of linear iterations in the algebraic solve, linear solver residual
   after each iteration etc.
   */
  std::unique_ptr<dealii::SolverControl> cn_;

  //! Hash table containing pointers of <a href="https://www.dealii.org/8.5.0/doxygen/deal.II/classdealii::PETScWrappers_1_1PreconditionBoomerAMG.html" style="color:blue"><b>BoomerAMG</b></a> preconditioner.
  std::unordered_map<int,
      std::unique_ptr<dealii::PETScWrappers::PreconditionBoomerAMG>> p_amg_;

  //! Hash table containing pointers of <a href="https://www.dealii.org/8.5.0/doxygen/deal.II/classdealii::PETScWrappers_1_1PreconditionBlockJacobi.html" style="color:blue"><b>block Jacobi</b></a> preconditioner.
  std::unordered_map<int,
      std::unique_ptr<dealii::PETScWrappers::PreconditionBlockJacobi>> p_bjacobi_;

  //! Hash table containing pointers of <a href="https://www.dealii.org/8.5.0/doxygen/deal.II/classdealii::PETScWrappers_1_1PreconditionParaSails.html" style="color:blue"><b>ParaSails</b></a> preconditioner.
  std::unordered_map<int,
      std::unique_ptr<dealii::PETScWrappers::PreconditionParaSails>> p_parasails_;

  //! Hash table containing pointers of <a href="https://www.dealii.org/8.5.0/doxygen/deal.II/classdealii::PETScWrappers_1_1PreconditionJacobi.html" style="color:blue"><b>Jacobi</b></a> preconditioner.
  std::unordered_map<int,
      std::unique_ptr<dealii::PETScWrappers::PreconditionJacobi>> p_jacobi_;

  //! Hash table containing pointers of <a href="https://www.dealii.org/8.5.0/doxygen/deal.II/classdealii::PETScWrappers_1_1PreconditionEisenstat.html" style="color:blue"><b>Eisenstat</b></a> (block SSOR) preconditioner.
  std::unordered_map<int,
      std::unique_ptr<dealii::PETScWrappers::PreconditionEisenstat>> p_eisenstat_;

  //! Hash table containing pointers of <a href="https://www.dealii.org/8.5.0/doxygen/deal.II/classdealii::PETScWrappers_1_1SparseDirectMUMPS.html" style="color:blue"><b>MUMPS</b></a> solver objects.
  std::unordered_map<int,
      std::unique_ptr<dealii::PETScWrappers::SparseDirectMUMPS>> direct_;
};

#endif //BART_SRC_EQUATION_LINEAR_ALGEBRA_H_
