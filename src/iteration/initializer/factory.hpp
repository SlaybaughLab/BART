#ifndef BART_SRC_ITERATION_INITIALIZER_FACTORY_HPP_
#define BART_SRC_ITERATION_INITIALIZER_FACTORY_HPP_

namespace bart::iteration::initializer {

enum class InitializerName {
  kInitializeFixedTermsOnce = 0,
  kInitializeFixedTermsAndResetMoments = 1,
};

} // namespace bart::iteration::initializer

#endif //BART_SRC_ITERATION_INITIALIZER_FACTORY_HPP_
