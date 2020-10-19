#include "instrumentation/instrument_array.hpp"

#include <deal.II/lac/vector.h>

#include "convergence/status.h"
#include "instrumentation/tests/instrument_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "utility/colors.h"

namespace  {

namespace convergence = bart::convergence;
namespace instrumentation = bart::instrumentation;
namespace utility = bart::utility;

template <typename InputType>
class InstrumentationInstrumentArrayTest : public ::testing::Test {
 public:
};

using InstrumentTypes = ::testing::Types<convergence::Status, std::pair<int, double>,
                                         std::pair<std::string, utility::Color>,
                                         dealii::Vector<double>>;

TYPED_TEST_SUITE(InstrumentationInstrumentArrayTest, InstrumentTypes);

} // namespace
