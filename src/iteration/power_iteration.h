#ifndef BART_SRC_ITERATION_POWER_ITERATION_H_
#define BART_SRC_ITERATION_POWER_ITERATION_H_

#include "eigen_base.h"

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
class PowerIteration : public EigenBase<dim> {
 public:
  /*!
   Class constructor.

   \param prm ParameterHandler object.
   */
  PowerIteration (const dealii::ParameterHandler &prm,
      std::shared_ptr<FundamentalData<dim>> &dat_ptr);

  //! Class destructor.
  ~PowerIteration ();

  /*!
   Eigenvalue iterations using power iteration scheme for eigenvalue calculations.

   \param equ_ptr Pointer of equation, i.e. EquationBase<dim> object.
   \return Void.
   */
  void EigenIterations (std::unique_ptr<EquationBase<dim>> &equ_ptr);
};

#endif //BART_SRC_ITERATION_POWER_ITERATION_H_
