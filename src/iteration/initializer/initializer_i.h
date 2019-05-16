#ifndef BART_SRC_ITERATION_INITIALIZER_INITIALIZER_I_H_
#define BART_SRC_ITERATION_INITIALIZER_INITIALIZER_I_H_

namespace bart {

namespace data {

struct System;

} // namespace data

namespace iteration {

namespace initializer {
/*! \brief Interface for class that initialize a system for an iteration.
 *
 */
class InitializerI {
 public:
  virtual ~InitializerI() = default;
  virtual void Initialize(data::System& sys) = 0;
};

} // namespace initializer

} // namespace iteration

} // namespace bart

#endif // BART_SRC_ITERATION_INITIALIZER_INITIALIZER_I_H_