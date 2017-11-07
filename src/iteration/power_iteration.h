#ifndef __power_iteration_h__
#define __power_iteration_h__

#include "eigen_base.h"

using namespace dealii;

//! This class provides power iteration scheme.
/*!
 This class implements power iteration for eigenvalue calculations. Assuming 
 \f$T\f$, \f$S\f$ and \f$F\f$ are transport, scattering and fission operators, 
 a power iteration scheme can be represented as
 \f[
 T\psi^{l+1}=S\psi^{l+1}+F\psi^l,
 \f]
 where \f$l\f$ stands for eigenvalue iteration index.
 
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
  
  /*!
   Eigenvalue iterations using power iteration scheme for eigenvalue calculations. 
   
   \param sflxes_proc Scalar fluxes for all groups living on current processor.
   \param equ_ptrs Pointers of equations, i.e. EquationBase<dim> objects.
   \param ig_ptr Pointer of in-group solver, i.e. IGBase<dim> object
   \param mg_ptr Pointer of MG solver, i.e. MGBase<dim> object.
   \return Void.
   */
  void eigen_iterations
  (std::vector<Vector<double> > &sflxes_proc,
   std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs,
   std_cxx11::shared_ptr<IGBase<dim> > ig_ptr,
   std_cxx11::shared_ptr<MGBase<dim> > mg_ptr);
};

#endif //__power_iteration_h__
