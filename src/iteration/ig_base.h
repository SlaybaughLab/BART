#ifndef __ig_base_h__
#define __ig_base_h__

#include "iteration_base.h"

using namespace dealii;

template <int dim>
class IGBase : public IterationBase<dim>
{
public:
  IGBase (const ParameterHandler &prm);
  virtual ~IGBase();
  
  // has to be provided
  virtual void solve_in_group
  (std::vector<Vector<double> > &sflxes_proc,
   std_cxx11::shared_ptr<EquationBase<dim> > equ_ptrs,
   unsigned int &g);
  
protected:
  const double err_phi_tol;
  
  Vector<double> sflx_proc_prev_ig;
};

template <int dim>
class SourceIteration : public IGBase<dim>
{
public:
  SourceIteration (const ParameterHandler &prm);
  ~SourceIteration ();
  
  void solve_in_group
  (std::vector<Vector<double> > &sflxes_proc,
   std_cxx11::shared_ptr<EquationBase<dim> > equ_ptr,
   unsigned int &g);
};

#endif //__ig_base_h__
