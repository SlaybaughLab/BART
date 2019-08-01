#ifndef BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_H_
#define BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_H_

#include "framework/builder/framework_builder_i.h"

namespace bart {

namespace framework {

namespace builder {

class FrameworkBuilder : public FrameworkBuilderI {
 public:
  FrameworkBuilder() = default;
  virtual ~FrameworkBuilder() = default;
};

} // namespace builder

} // namespace framework

} // namespace bart

#endif //BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_H_
