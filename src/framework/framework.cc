#include "framework/framework.h"

namespace bart {

namespace framework {

Framework::Framework(std::unique_ptr<system::System> system_ptr)
    : system_ptr_(std::move(system_ptr)) {}

} // namespace framework

} // namespace bart

