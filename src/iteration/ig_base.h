#ifndef BART_SRC_ITERATION_IG_BASE_H_
#define BART_SRC_ITERATION_IG_BASE_H_

#include "iteration_base.h"

using namespace dealii;

//! This class provides in-group solving scheme.
/*!
 This class serves as the base class for in group solvers. The main (only)
 functionality is to provide an abstract function for solve in a specific group.
 Overriding has to be provided in derived class.

 \author Weixiong Zheng
 \date 2017/09
 */
template <int dim>
class IGBase : public IterationBase<dim> {
 public:
  /*!
   Class constructor.

   \param prm Const ParameterHandler object.
   */
  IGBase (const ParameterHandler &prm,
      std::shared_ptr<FundamentalData<dim>> dat_ptr);

  //! Class destructor.
  virtual ~IGBase();

  /*!
   Abstract function for in group solving functionality. For instance, one could
   provid overriding as one-pass solve for diffusion/NDA, or iteration based in
   group solving such as source iteration.

   \param sflxes_proc A vector of all scalar fluxes on current processor.
   \param equ_ptr A pointer to the equation of interest.
   \param g Group index.
   */
  virtual void IGIterations (
      std::uniq_ptr<EquationBase<dim>> equ_ptr,
      const int &g);

 protected:
  const double err_phi_tol_;//!< Tolerance used in iterative in group scheme for convergence check.
};

#endif //BART_SRC_ITERATION_IG_BASE_H_
