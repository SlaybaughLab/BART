#ifndef __power_iteration_h__
#define __power_iteration_h__

#include "eigen_base.h"

template <int dim>
class PowerIteration : public EigenBase<dim>
{
public:
  PowerIteration (ParameterHandler);
  ~PowerIteration ();
  
  void eigen_iterations
  (std::vector<Vector<double> > &sflxes_proc,
   std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > equ_ptrs);
}

#endif //__power_iteration_h__
