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
  PreconditionerSolver (ParameterHandler &prm);
  ~PreconditionerSolver ();

private:

};

#endif //__preconditioner_solver_h__
