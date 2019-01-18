#ifndef BART_SRC_PROBLEM_LOCATOR_H_
#define BART_SRC_PROBLEM_LOCATOR_H_

#include <memory>
#include "../problem/parameters_i.h"

namespace bart {

namespace problem {

class Locator {
 public:
  static ParametersI *GetParameters() { return parameters_; }

  static void Provide(ParametersI *parameters)
  {
    parameters_ = parameters;
  }
  
 private:
  static ParametersI *parameters_;
};

} // namespace problem

} // namespace bart

#endif // BART_SRC_PROBLEM_LOCATOR_H_
