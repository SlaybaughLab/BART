#ifndef BART_SRC_ITERATION_SOURCE_ITERATION_H_
#define BART_SRC_ITERATION_SOURCE_ITERATION_H_

#include "ig_base.h"

//! This class provides source iteration scheme for in group solve.
/*!
 \author Weixiong Zheng
 \date 2017/08
 */
template <int dim>
class SourceIteration : public IGBase<dim> {
 public:
  /*!
   Class constructor.

   \param prm Const ParameterHandler object.
   */
  SourceIteration (const dealii::ParameterHandler &prm,
      std::shared_ptr<FundamentalData<dim>> &dat_ptr);

  //! Class destructor.
  ~SourceIteration ();

  /*!
   A function to solve SN in group using source iteration scheme.

   \param sflxes_proc Scalar fluxes for all groups living on current processor.
   \param equ_ptr A pointer to the equation of interest.
   \param g Group index.
   */
  void IGIterations ( std::unique_ptr<EquationBase<dim>> &equ_ptr,
      const int &g);
};

#endif //BART_SRC_ITERATION_SOURCE_ITERATION_H_
