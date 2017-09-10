#ifndef __mg_base_h__
#define __mg_base_h__

#include "../transport/equation_base.h"
#include "in_group_base.h"

using namespace dealii;

template <int dim>
class MGBase : public IterationBase<dim>
{
public:
  MGBase (ParameterHandler &prm);
  virtual ~MGBase ();
  
  virtual void do_iterations
  (std::vector<Vector<double> > &sflxes_proc,
   std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs);
  
  virtual void mg_iterations
  (std::vector<Vector<double> > &sflxes_proc,
   std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs);
protected:
  const double err_phi_tol;
  
  std_cxx11::shared_ptr<InGroupBase<dim> > ig_ptr;
  
  std::vector<Vector<double> > sflxes_proc_prev_mg;
}

#endif//__mg_base_h__
