#ifndef __gauss_seidel_h__
#define __gauss_seidel_h__

#include "mg_base.h"

//! This class provides Gauss-Seidel scheme for MG calculations.
/*!
 This class implements Gauss-Seidel scheme for MG calculations. Specifically,
 MGBase<dim>::nonthermal_solves and MGBase<dim>::thermal_iterations are overriden
 using Gauss-Seidel iterations.

 A Gauss-Seidel iteration can be represented as:
 \f[
 T_g\psi_g^{l+1}=\sum\limits_{g'=0}^{g}S_{g'\to g}\psi_{g'}^{l+1}+\sum\limits_{g'=g+1}^{G-1}S_{g'\to g}\psi_{g'}^{l}+Q,
 \f]
 where \f$T\f$, \f$S\f$ and \f$Q\f$ are transport operator, scattering operator
 and source and \f$l\f$ stands for MG iteration index.

 \author Weixiong Zheng
 \date 2017/09
 */
template <int dim>
class GaussSeidel : public MGBase<dim> {
 public:
  /*!
   Class constructor.

   \param prm const ParameterHandler object.
   */
  GaussSeidel (const ParameterHandler &prm,
      std::shared_ptr<FundamentalData<dim>> dat_ptr);

  //! Class destructor.
  ~GaussSeidel ();

  /*!
   A function overriding MGBase<dim>::nonthermal_solves. The function implements
   Gauss-Seidel scheme with one-pass serial solves (energy groupwise) from highest
   energy group until reaching out to the last group before the thermal group
   boundary.

   \param equ_ptrs A vector of pointers of EquationBase<dim> object.
   \return Void.
   */
  void NonthermalSolves (std::unique_ptr<EquationBase<dim>> equ_ptr);

  /*!
   A function overriding MGBase<dim>::thermal_iterations. The function implements
   Gauss-Seidel scheme for solving transport problems in thermal groups
   iteratively. The starting group is the highest energy group affected by
   upscattering.

   \param equ_ptrs A vector of pointers of EquationBase<dim> object.
   \return Void.
   */
  void ThermalIterations (std::unique_ptr<EquationBase<dim>> equ_ptr);
};

#endif //__gauss_Seidel_h__
