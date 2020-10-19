#include "instrumentation/instrument_array.hpp"

#include <deal.II/lac/vector.h>

#include "convergence/status.h"
#include "instrumentation/tests/instrument_mock.h"
#include "system/moments/spherical_harmonic_i.h"
#include "test_helpers/gmock_wrapper.h"
#include "utility/colors.h"

namespace  {

namespace convergence = bart::convergence;
namespace instrumentation = bart::instrumentation;
namespace utility = bart::utility;
namespace moments = bart::system::moments;


template <typename InputType>
class InstrumentationInstrumentArrayTest : public ::testing::Test {
 public:
};

using InstrumentTypes = ::testing::Types<convergence::Status, std::pair<int, double>,
                                         std::pair<std::string, utility::Color>,
                                         dealii::Vector<double>,
                                         moments::SphericalHarmonicI>;

TYPED_TEST_SUITE(InstrumentationInstrumentArrayTest, InstrumentTypes);

TYPED_TEST(InstrumentationInstrumentArrayTest, Dummy) {
  EXPECT_TRUE(true);
}



} // namespace
