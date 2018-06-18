#ifndef BART_SRC_ITERATION_MG_BASE_H_
#define BART_SRC_ITERATION_MG_BASE_H_

#include "iteration_base.h"
#include "ig_base.h"

//! This class provides abstract functionalities for MG calculations.
/*!
 This class serves as the base class of MG iterations. It inherits from
 IterationBase.

 \author Weixiong Zheng
 \date 2017/08
 */
template <int dim>
class MGBase : public IterationBase<dim> {
 public:
  /*!
   A constructor of MGBase.

   \param prm A ParameterHandler object containing all the parameters needed.
   */
  MGBase (const ParameterHandler &prm,
      std::shared_ptr<FundamentalData<dim>> dat_ptr);

  //! A virtual destructor for MGBase
  virtual ~MGBase ();

  /*!
   This function will call mg_iterations to perform MG iterations in fixed source
   problems. This function shall not be called in eigenvalue calculations.

   \todo Implement NDA scheme in this function.

   \param sflxes_proc A vector of scalar fluxes for all groups.
   \param equ_ptrs A vector of shared_ptr's of EquationBase objects.
   \param ig_ptr A shared_ptr of IGBase<dim> object.
   \return Void.
   */
  virtual void DoIterations (std::unordered_map<std::string,
      std::unique_ptr<EquationBase<dim>>> &equ_ptrs);

  /*!
   A virtual function performing MG iterations. By default, it will call
   nonthermal_solves for nonthermal MG calculations and thereafter the
   thermal_iterations to iterate on the upscattering.

   Overriding could be provided for this virtual function if nonthermal and
   thermal solves are not separated. An example is Jacobian Free Newton Krylov
   method.

   There are two case of calling this function: one is to call this function inside
   eigen iterations and the other is in do_iterations, which is designed for MG
   fixed source problem calculations.

   \param sflxes_proc A vector of scalar fluxes for all groups living on current
   processor.
   \param equ_ptrs A vector of shared_ptr's of EquationBase objects.
   \param ig_ptr A shared_ptr of IGBase object.
   \return Void.
   */
  virtual void MGIterations (std::unique_ptr<EquationBase<dim>> equ_ptr);

  /*!
   This virtual function performs energy solves over nonthermal groups.

   Usually, nonthermal groups have no upscattering. So this function is a group-by-
   group one-pass solving until reaching the thermal group. It will not be called
   if algorithms like JFNK are called

   \param equ_ptr Pointer of EquationBase object.
   \return Void.
   */
  virtual void NonthermalSolves(std::unique_ptr<EquationBase<dim>> equ_ptr);

  /*!
   This virtual function performs iterative energy solves over thermal groups.

   Thermal groups have upscattering for applications like LWR. So this function is
   to solve for thermal groups iteratively. It will not be called if algorithms like
   JFNK are called.

   \param equ_ptrs A vector of shared_ptr's of EquationBase objects.
   \return Void.
   */
  virtual void ThermalIterations (std::unique_ptr<EquationBase<dim>> equ_ptr);

 protected:
  unsigned int g_thermal_;//!< Starting group index where upscattering exists.
  const double err_phi_tol_;//!< Multigroup iteration tolerance for convergence check.
};

#endif //BART_SRC_ITERATION_MG_BASE_H_
