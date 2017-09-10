#ifndef __in_group_base_h__
#define __in_group_base_h__

#include "in_group_base.h"

template <int dim>
class InGroupBase : public IterationBase<dim>
{
public:
  InGroupBase (ParameterHandler &prm);
  virtual ~InGroupBase();
  
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
class SourceIteration : public InGroupBase<dim>
{
public:
  SourceIteration (ParameterHandler &prm);
  ~SourceIteration ();
  
  void solve_in_group
  (std::vector<Vector<double> > &sflxes_proc,
   std_cxx11::shared_ptr<EquationBase<dim> > equ_ptrs,
   unsigned int &g);
};

#endif //__in_group_base_h__
