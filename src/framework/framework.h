#ifndef BART_SRC_FRAMEWORK_FRAMEWORK_H_
#define BART_SRC_FRAMEWORK_FRAMEWORK_H_

#include "framework/framework_i.h"

namespace framework {

class Framework : public FrameworkI {
 public:
  virtual ~Framework() = default;
};

} // namespace framework

#endif //BART_SRC_FRAMEWORK_FRAMEWORK_H_
