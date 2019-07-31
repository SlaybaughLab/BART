#include "framework/framework.h"

namespace bart {

namespace framework {

Framework::Framework(
    std::unique_ptr<system::System> system_ptr,
    std::unique_ptr<Initializer> initializer_ptr,
    std::unique_ptr<OuterIterator> outer_iterator_ptr)
    : system_ptr_(std::move(system_ptr)),
      initializer_ptr_(std::move(initializer_ptr)),
      outer_iterator_ptr_(std::move(outer_iterator_ptr))
    {}

} // namespace framework

} // namespace bart

