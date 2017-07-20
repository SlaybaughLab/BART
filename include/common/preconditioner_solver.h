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
                        unsigned int &n_total_ho_vars,
                        MPI_Comm &mpi_communicator);
  ~PreconditionerSolver ();
  
  // HO solver related member functions
  void initialize_ho_preconditioners
  (std::vector<PETScWrappers::MPI::SparseMatrix*> &ho_syses,
   std::vector<PETScWrappers::MPI::Vector*> &ho_rhses);
  
  void ho_solve (std::vector<PETScWrappers::MPI::SparseMatrix*> &ho_syses,
                 std::vector<PETScWrappers::MPI::Vector*> &ho_psis,
                 std::vector<PETScWrappers::MPI::Vector*> &ho_rhses);
  
  // NDA solver related member functions
  void reinit_nda_preconditioners
  (std::vector<PETScWrappers::MPI::SparseMatrix*> &ho_syses,
   std::vector<PETScWrappers::MPI::Vector*> &ho_rhses);
  
  void nda_solve
  (PETScWrappers::MPI::SparseMatrix &nda_sys,
   PETScWrappers::MPI::Vector &nda_phi,
   PETScWrappers::MPI::Vector &nda_rhs,
   unsigned int &g);
  
private:
  const unsigned int n_group;
  const unsigned int n_total_ho_vars;
  const bool do_nda;
  
  bool have_reflective_bc;
  double ho_ssor_omega;
  double nda_ssor_omega;
  
  std::string transport_model_name;
  std::string ho_linear_solver_name;
  std::string ho_preconditioner_name;
  std::string nda_linear_solver_name;
  std::string nda_preconditioner_name;
  
  std::vector<bool> ho_direct_init;
  std::vector<bool> nda_direct_init;
  std::vector<unsigned int> ho_linear_iters;
  std::vector<unsigned int> nda_linear_iters;
  
  // HO solver related variables
  std::vector<std_cxx11::shared_ptr<SolverControl> > ho_cn;
  std::vector<std_cxx11::shared_ptr<PETScWrappers::PreconditionBoomerAMG> > pre_ho_amg;
  std::vector<std_cxx11::shared_ptr<PETScWrappers::PreconditionBlockJacobi> > pre_ho_bjacobi;
  std::vector<std_cxx11::shared_ptr<PETScWrappers::PreconditionParaSails> > pre_ho_parasails;
  std::vector<std_cxx11::shared_ptr<PETScWrappers::PreconditionJacobi> > pre_ho_jacobi;
  std::vector<std_cxx11::shared_ptr<PETScWrappers::PreconditionEisenstat> > pre_ho_eisenstat;
  std::vector<std_cxx11::shared_ptr<PETScWrappers::SparseDirectMUMPS> > ho_direct;
  
  // NDA solver related variables
  std::vector<std_cxx11::shared_ptr<SolverControl> > nda_cn;
  std::vector<std_cxx11::shared_ptr<PETScWrappers::PreconditionBoomerAMG> > pre_nda_amg;
  std::vector<std_cxx11::shared_ptr<PETScWrappers::PreconditionBlockJacobi> > pre_nda_bjacobi;
  std::vector<std_cxx11::shared_ptr<PETScWrappers::PreconditionParaSails> > pre_nda_parasails;
  std::vector<std_cxx11::shared_ptr<PETScWrappers::PreconditionJacobi> > pre_nda_jacobi;
  std::vector<std_cxx11::shared_ptr<PETScWrappers::PreconditionEisenstat> > pre_nda_eisenstat;
  std::vector<std_cxx11::shared_ptr<PETScWrappers::SparseDirectMUMPS> > nda_direct;
  
  MPI_Comm mpi_communicator;
};

#endif //__preconditioner_solver_h__
