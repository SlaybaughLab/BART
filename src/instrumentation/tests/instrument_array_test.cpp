#include "instrumentation/instrument_array.hpp"

#include <deal.II/lac/vector.h>

#include "convergence/status.h"
#include "instrumentation/tests/instrument_mock.h"
#include "system/moments/spherical_harmonic_i.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"
#include "utility/colors.h"

namespace  {

namespace convergence = bart::convergence;
namespace instrumentation = bart::instrumentation;
namespace moments = bart::system::moments;
namespace utility = bart::utility;
namespace test_helpers = bart::test_helpers;

template <typename InputType>
class InstrumentationInstrumentArrayTest : public ::testing::Test {
 public:
  using InstrumentType = instrumentation::InstrumentMock<InputType>;
  instrumentation::InstrumentArray<InputType> test_instrument_array_;
};

using InstrumentTypes = ::testing::Types<moments::SphericalHarmonicI>;

TYPED_TEST_SUITE(InstrumentationInstrumentArrayTest, InstrumentTypes);

TYPED_TEST(InstrumentationInstrumentArrayTest, AddInstrument) {
  using InstrumentType = typename InstrumentationInstrumentArrayTest<TypeParam>::InstrumentType;
  auto& test_instrument_array = this->test_instrument_array_;

  const int n_instruments{ test_helpers::RandomInt(2, 5) };
  EXPECT_EQ(test_instrument_array.size(), 0);
  for (int i = 0; i < n_instruments; ++i) {
    EXPECT_NO_THROW({
      test_instrument_array.AddInstrument(std::make_unique<InstrumentType>());
    });
  }
  EXPECT_EQ(test_instrument_array.size(), n_instruments);
}



} // namespace
