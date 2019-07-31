#ifndef BART_SRC_FRAMEWORK_FRAMEWORK_H_
#define BART_SRC_FRAMEWORK_FRAMEWORK_H_

#include <memory>

#include "system/system.h"
#include "framework/framework_i.h"

namespace bart {

namespace framework {

class Framework : public FrameworkI {
 public:
  Framework(std::unique_ptr<system::System> system_ptr);
  virtual ~Framework() = default;

  system::System* system() const {
    return system_ptr_.get();
  }

 protected:
  std::unique_ptr<system::System> system_ptr_ = nullptr;

};

} // namespace framework

} // namespace bart

#endif //BART_SRC_FRAMEWORK_FRAMEWORK_H_
