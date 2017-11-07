#ifndef __power_iteration_h__
#define __power_iteration_h__

#include "eigen_base.h"

using namespace dealii;

//! This class provides power iteration scheme.
/*!
 This class implements power iteration for eigenvalue calculations.
 
 \author Weixiong Zheng
 \date 2017/09
 */
template <int dim>
class PowerIteration : public EigenBase<dim>
{
public:
  /*!
   Class constructor.
   
   \param prm ParameterHandler object.
   */
  PowerIteration (const ParameterHandler &prm);
  
  //! Class destructor.
  ~PowerIteration ();
  
  void eigen_iterations
  (std::vector<Vector<double> > &sflxes_proc,
   std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs,
   std_cxx11::shared_ptr<IGBase<dim> > ig_ptr,
   std_cxx11::shared_ptr<MGBase<dim> > mg_ptr);
};

#endif //__power_iteration_h__
