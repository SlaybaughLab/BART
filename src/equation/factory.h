#ifndef BART_SRC_EQUATION_FACTORY_H_
#define BART_SRC_EQUATION_FACTORY_H_

#include <memory>

#include "equation_base.h"
#include "../common/computing_data.h"

namespace bart {

namespace equation {

template <int dim>
class Factory {
 private:
  Factory() = default;
  ~Factory() = default;
 public:
  static std::unique_ptr<EquationBase<dim>> CreateEquation(
      const dealii::ParameterHandler &prm,
      std::shared_ptr<FundamentalData<dim>> &dat_ptr);

};
  
} // namespace equation

} // namespace bart

#endif // BART_SRC_EQUATION_FACTORY_H_
