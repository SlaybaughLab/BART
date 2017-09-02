#ifndef __mg_base_h__
#define __mg_base_h__

#include "../transport/equation_base.h"

using namespace dealii;

template <int dim>
class MGBase : public IterationBase<dim>
{
public:
  MGBase ();
  virtual ~MGBase ();
  
  virtual void do_iterations
  (std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs);
  
  virtual void mg_iterations ();
  
  virtual void generate_system_matrices
  (std::vector<PETScWrappers::MPI::SparseMatrix*> &sys_mats);
  
  virtual void generate_group_rhses
  (std::vector<PETScWrappers::MPI::Vector*> &group_rhses, unsigned int &g);
  
  virtual void iteration_over_groups
  (std::vector<PETScWrappers::MPI::Vector*> &group_rhses);
  
protected:
  const double err_phi_tol;
  
  std::vector<Vector<double> > sflx_proc_old;
  // auxiliary flux if it's HOLO methodology: it will be for HO flux
  std::vector<Vector<double> > sflx_proc_aux;
}

#endif//__mg_base_h__
