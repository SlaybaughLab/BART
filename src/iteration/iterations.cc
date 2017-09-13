#include <deal.II/base/index_set.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/cell_id.h>

#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/solver_bicgstab.h>

#include <algorithm>

#include "iterations.h"
#include "../common/bart_tools.h"

using namespace dealii;

template <int dim>
Iterations<dim>::Iterations
(const ParameterHandler &prm,
 const DoFHandler<dim> &dof_handler
 const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
 const std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr,
 const std_cxx11::shared_ptr<MaterialProperties> mat_ptr)
:
transport_name(prm.get("transport model")),
is_eigen_problem(prm.get_bool("do eigenvalue calculations"))
{
  mat_ptr = build_material (prm);
  // vectors containing all the pointers to equations. Size is 2 if NDA is used.
  bool do_nda = prm.get_bool ("do NDA");
  equ_ptrs.resize (do_nda?2:1);
  equ_ptrs[0] = build_equation (transport_name, prm, dof_handler, msh_ptr, aqd_ptr, mat_ptr);
  if (do_nda)
    equ_ptrs[1] = build_equation ("nda", prm, dof_handler, msh_ptr, aqd_ptr, mat_ptr);
}

template <int dim>
Iterations<dim>::~Iterations ()
{
}

template <int dim>
void Iterations<dim>::initialize_assembly_related_objects
(FE_Poly<TensorProductPolynomials<dim>,dim,dim>* fe)
{
  for (unsigned int i=0; i<equ_ptrs.size(); ++i)
    equ_ptrs[i]->initialize_assembly_related_objects (fe);
}

template <int dim>
void Iterations<dim>::initialize_system_matrices_vectors
(SparsityPatternType &dsp,
 IndexSet &local_dofs,
 std::vector<Vector<double> > &sflxes_proc)
{
  // 1. each equation pointer contains system matrices and vectors from PETSc, per se
  // This function is to tell the shapes and sizes of those vectors and matrices
  
  // 2. initialize sflxes on this proc with unit values in the right shapes
  for (unsigned int i=0; i<equ_ptrs.size(); ++i)
    equ_ptrs[i]->initialize_system_matrices_vectors (dsp, local_dofs, sflxes_proc);
}

template <int dim>
void Iterations<dim>::solve_problems (std::vector<Vector<double> > &sflxes_proc)
{
  if (is_eigen_problem)
  {
    std_cxx11::shared_ptr<EigenBase<dim> > pro_ptr = build_eigen_problem (prm);
    pro_ptr->do_iterations (sflxes_proc, equ_ptrs);
    pro_ptr->get_keff (keff);
  }
  else
  {
    std_cxx11::shared_ptr<MGBase<dim> > pro_ptr = build_mg_problem (prm);
    pro_ptr->do_iterations (sflxes_proc, equ_ptrs);
  }
}

template <int dim>
void Iterations<dim>::get_keff (double &k)
{
  AssertThrow (is_eigen_problem,
               ExcMessage("Only eigen problems have keff"));
  k = keff;
}

// explicit instantiation to avoid linking error
template class Iterations<2>;
template class Iterations<3>;
