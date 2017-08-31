#ifndef __preconditioner_solver_h__
#define __preconditioner_solver_h__

#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/base/parameter_handler.h>

#include <vector>
#include <string>

using namespace dealii;

class PreconditionerSolver
{
public:
  PreconditionerSolver (ParameterHandler &prm,
                        std::string equation_name,
                        unsigned int& n_total_vars);
  ~PreconditionerSolver ();
  
  // HO solver related member functions
  void initialize_preconditioners
  (std::vector<PETScWrappers::MPI::SparseMatrix*> &sys_mats,
   std::vector<PETScWrappers::MPI::Vector*> &sys_rhses);
  
  void linear_algebra_solve
  (std::vector<PETScWrappers::MPI::SparseMatrix*> &sys_mats,
   std::vector<PETScWrappers::MPI::Vector*> &sys_flxes,
   std::vector<PETScWrappers::MPI::Vector*> &sys_rhses,
   unsigned int &i/*component number*/);

private:
  const unsigned int n_total_vars;
  const std::string equation_name;
  
  bool have_reflective_bc;
  double ho_ssor_omega;
  double nda_ssor_omega;
  
  std::string linear_solver_name;
  std::string preconditioner_name;
  
  std::vector<bool> ho_direct_init;
  std::vector<bool> nda_direct_init;
  std::vector<unsigned int> ho_linear_iters;
  std::vector<unsigned int> nda_linear_iters;
  
  // solver related variables
  // TODO: Add some other PETSc preconditioners existing in deal.II
  // solver control pointers: handling information for/from linear solvers
  std::vector<std_cxx11::shared_ptr<SolverControl> > cn;
  // preconditioner pointers: different preconditioners for iterative linear solvers
  std::vector<std_cxx11::shared_ptr<PETScWrappers::PreconditionBoomerAMG> > pre_amg;
  std::vector<std_cxx11::shared_ptr<PETScWrappers::PreconditionBlockJacobi> > pre_bjacobi;
  std::vector<std_cxx11::shared_ptr<PETScWrappers::PreconditionParaSails> > pre_parasails;
  std::vector<std_cxx11::shared_ptr<PETScWrappers::PreconditionJacobi> > pre_jacobi;
  std::vector<std_cxx11::shared_ptr<PETScWrappers::PreconditionEisenstat> > pre_eisenstat;
  // direct solver pointers: in case iterative solvers are not chosen, one my use
  // direct solver based MUMPS
  std::vector<std_cxx11::shared_ptr<PETScWrappers::SparseDirectMUMPS> > direct;
};

#endif //__preconditioner_solver_h__
