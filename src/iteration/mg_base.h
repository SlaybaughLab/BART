#ifndef __mg_base_h__
#define __mg_base_h__

#include "iteration_base.h"
#include "ig_base.h"

using namespace dealii;

template <int dim>
class MGBase : public IterationBase<dim>
{
public:
  MGBase (const ParameterHandler &prm);
  virtual ~MGBase ();
  
  virtual void do_iterations
  (std::vector<Vector<double> > &sflxes_proc,
   std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs,
   std_cxx11::shared_ptr<IGBase<dim> > ig_ptr);
  
  virtual void mg_iterations
  (std::vector<Vector<double> > &sflxes_proc,
   std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs,
   std_cxx11::shared_ptr<IGBase<dim> > ig_ptr);
  
  virtual void nonthermal_solves
  (std::vector<Vector<double> > &sflxes_proc,
   std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs,
   std_cxx11::shared_ptr<IGBase<dim> > ig_ptr);
  
  virtual void thermal_iterations
  (std::vector<Vector<double> > &sflxes_proc,
   std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs,
   std_cxx11::shared_ptr<IGBase<dim> > ig_ptr);
protected:
  const double err_phi_tol;
  
  unsigned int g_thermal;
  
  std::vector<Vector<double> > sflxes_proc_prev_mg;
};

#endif//__mg_base_h__
