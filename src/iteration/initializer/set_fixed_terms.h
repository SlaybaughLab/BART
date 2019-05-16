#ifndef BART_SRC_ITERATION_INITIALIZER_SET_FIXED_TERMS_H_
#define BART_SRC_ITERATION_INITIALIZER_SET_FIXED_TERMS_H_

#include "iteration/initializer/initializer_i.h"

namespace bart {

namespace iteration {

namespace initializer {
/*! \brief Initializes a system by setting the fixed terms.
 *
 * This class takes a class of type iteration::updater::FixedUpdaterI and
 * uses the UpdateFixedTerms method to initialize the system.
 *
 */
class SetFixedTerms : public InitializerI {
 public:
  virtual ~SetFixedTerms() = default;
};

} // namespace initializer

} // namespace iteration

} // namespace bart

#endif // BART_SRC_ITERATION_INITIALIZER_SET_FIXED_TERMS_H_