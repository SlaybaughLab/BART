#ifndef __power_iteration_h__
#define __power_iteration_h__

#include "eigen_base.h"

using namespace dealii;

template <int dim>
class PowerIteration : public EigenBase<dim>
{
public:
  PowerIteration (const ParameterHandler &prm);
  ~PowerIteration ();
  
  void do_iterations
  (std::vector<Vector<double> > &sflxes_proc,
   std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs,
   std_cxx11::shared_ptr<IGBase<dim> > ig_ptr,
   std_cxx11::shared_ptr<MGBase<dim> > mg_ptr);
  
  void eigen_iterations
  (std::vector<Vector<double> > &sflxes_proc,
   std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs,
   std_cxx11::shared_ptr<IGBase<dim> > ig_ptr,
   std_cxx11::shared_ptr<MGBase<dim> > mg_ptr);
};

#endif //__power_iteration_h__
