#ifndef BART_SRC_FRAMEWORK_FRAMEWORK_I_H_
#define BART_SRC_FRAMEWORK_FRAMEWORK_I_H_

namespace bart {

namespace framework {

class FrameworkI {
 public:
  virtual ~FrameworkI() = default;
  virtual void SolveSystem() = 0;
};

} // namespace framework

} // namespace bart

#endif //BART_SRC_FRAMEWORK_FRAMEWORK_I_H_
