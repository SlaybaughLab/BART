#ifndef BART_SRC_FRAMEWORK_FRAMEWORK_HPP_
#define BART_SRC_FRAMEWORK_FRAMEWORK_HPP_

#include <memory>
#include <ostream>

#include "iteration/outer/outer_iteration_i.hpp"
#include "iteration/initializer/initializer_i.h"
#include "results/output_i.h"
#include "system/system.h"
#include "framework/framework_i.h"

namespace bart {

namespace framework {

class Framework : public FrameworkI {
 public:
  using Initializer = iteration::initializer::InitializerI;
  using OuterIterator = iteration::outer::OuterIterationI;
  using ResultsOutput = results::OutputI;

  Framework(
      std::unique_ptr<system::System> system_ptr,
      std::unique_ptr<Initializer> initializer_ptr,
      std::unique_ptr<OuterIterator> outer_iterator_ptr,
      std::unique_ptr<ResultsOutput> results_output_ptr = nullptr);
  virtual ~Framework() = default;

  void SolveSystem() override;

  void OutputResults(std::ostream& output_stream);

  void OutputMasterFile(std::ostream& output_stream,
                        const std::vector<std::string>& filenames,
                        const int process_id);

  Initializer* initializer_ptr() const {
    return initializer_ptr_.get();
  }

  system::System* system() const override {
    return system_ptr_.get();
  }

  OuterIterator* outer_iterator_ptr() const {
    return outer_iterator_ptr_.get();
  }

  ResultsOutput* results_output_ptr() const {
    return results_output_ptr_.get();
  }

 protected:
  std::unique_ptr<system::System> system_ptr_ = nullptr;
  std::unique_ptr<Initializer> initializer_ptr_ = nullptr;
  std::unique_ptr<OuterIterator> outer_iterator_ptr_ = nullptr;
  std::unique_ptr<ResultsOutput> results_output_ptr_ = nullptr;
};

} // namespace framework

} // namespace bart

#endif //BART_SRC_FRAMEWORK_FRAMEWORK_HPP_
