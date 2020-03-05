#ifndef BART_SRC_ITERATION_INITIALIZER_INITIALIZE_FIXED_TERMS_H_
#define BART_SRC_ITERATION_INITIALIZER_INITIALIZE_FIXED_TERMS_H_

#include "iteration/initializer/initializer_i.h"

namespace bart {

namespace iteration {

namespace initializer {

class InitializeFixedTerms : public InitializerI {
 public:
  void Initialize(system::System &sys) override;
  virtual ~InitializeFixedTerms() = default;
};

} // namespace initializer

} // namespace iteration

} // namespace bart

#endif //BART_SRC_ITERATION_INITIALIZER_INITIALIZE_FIXED_TERMS_H_
