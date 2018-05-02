#ifndef BART_SRC_EQUATION_EQUATION_BASE_H_
#define BART_SRC_EQUATION_EQUATION_BASE_H_

#include "../common/fundamental_data.h"

template <int dim>
class EquationBase {
 public:
  EquationBase (dealii::ParameterHandler &prm,
      std::unique_ptr<FundamentalData<dim>> &fund_dat_ptr);
  ~EquationBase ();

 protected:
  std::unique_ptr<FundamentalData<dim>> dat_ptr;
};

#endif //BART_SRC_EQUATION_EQUATION_BASE_H_
