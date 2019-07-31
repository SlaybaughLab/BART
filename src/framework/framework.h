#ifndef BART_SRC_FRAMEWORK_FRAMEWORK_H_
#define BART_SRC_FRAMEWORK_FRAMEWORK_H_

#include <memory>

#include "iteration/initializer/initializer_i.h"
#include "system/system.h"
#include "framework/framework_i.h"

namespace bart {

namespace framework {

class Framework : public FrameworkI {
 public:
  using Initializer = iteration::initializer::InitializerI;

  Framework(
      std::unique_ptr<system::System> system_ptr,
      std::unique_ptr<Initializer> initializer_ptr);
  virtual ~Framework() = default;

  Initializer* initializer_ptr() const {
    return initializer_ptr_.get();
  }

  system::System* system() const {
    return system_ptr_.get();
  }

 protected:
  std::unique_ptr<system::System> system_ptr_ = nullptr;
  std::unique_ptr<Initializer> initializer_ptr_ = nullptr;

};

} // namespace framework

} // namespace bart

#endif //BART_SRC_FRAMEWORK_FRAMEWORK_H_
