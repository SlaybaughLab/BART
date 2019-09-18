#include "framework/framework.h"

namespace bart {

namespace framework {

Framework::Framework(
    std::unique_ptr<system::System> system_ptr,
    std::unique_ptr<Initializer> initializer_ptr,
    std::unique_ptr<OuterIterator> outer_iterator_ptr,
    std::unique_ptr<ResultsOutput> results_output_ptr)
    : system_ptr_(std::move(system_ptr)),
      initializer_ptr_(std::move(initializer_ptr)),
      outer_iterator_ptr_(std::move(outer_iterator_ptr)),
      results_output_ptr_(std::move(results_output_ptr)) {

  AssertThrow(system_ptr_ != nullptr,
              dealii::ExcMessage("System pointer passed to "
                                 "Framework constructor is null"));
  AssertThrow(initializer_ptr_!= nullptr,
              dealii::ExcMessage("Initializer pointer passed to "
                                 "Framework constructor is null"));
  AssertThrow(outer_iterator_ptr_ != nullptr,
              dealii::ExcMessage("Outer iterator pointer passed to "
                                 "Framework constructor is null"));

}

void Framework::SolveSystem() {
  initializer_ptr_->Initialize(*system_ptr_);
  outer_iterator_ptr_->IterateToConvergence(*system_ptr_);
}

} // namespace framework

} // namespace bart

