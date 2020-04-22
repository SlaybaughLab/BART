#ifndef BART_SRC_ITERATION_INITIALIZER_INITIALIZER_I_H_
#define BART_SRC_ITERATION_INITIALIZER_INITIALIZER_I_H_

#include "system/system.h"
#include "utility/has_description.h"

namespace bart {

namespace iteration {

namespace initializer {
/*! \brief Interface for class that initialize a system for an iteration.
 *
 */
class InitializerI : public utility::HasDescription {
 public:
  virtual ~InitializerI() = default;
  virtual void Initialize(system::System& sys) = 0;
};

} // namespace initializer

} // namespace iteration

} // namespace bart

#endif // BART_SRC_ITERATION_INITIALIZER_INITIALIZER_I_H_