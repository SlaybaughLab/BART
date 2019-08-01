#ifndef BART_SRC_FRAMEWORK_BUILDER_CFEM_FRAMEWORK_BUILDER_H_
#define BART_SRC_FRAMEWORK_BUILDER_CFEM_FRAMEWORK_BUILDER_H_

#include "framework/builder/framework_builder_i.h"

namespace bart {

namespace framework {

namespace builder {

class CFEM_FrameworkBuilder : public FrameworkBuilderI {
 public:
  CFEM_FrameworkBuilder() = default;
  virtual ~CFEM_FrameworkBuilder() = default;
};

} // namespace builder

} // namespace framework

} // namespace bart

#endif //BART_SRC_FRAMEWORK_BUILDER_CFEM_FRAMEWORK_BUILDER_H_
