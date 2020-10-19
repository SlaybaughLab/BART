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

using ::testing::WhenDynamicCastTo, ::testing::NotNull;

template <typename InputType>
class InstrumentationInstrumentArrayTest : public ::testing::Test {
 public:
  using InstrumentType = instrumentation::InstrumentMock<InputType>;
  using InstrumentArrayType = instrumentation::InstrumentArray<InputType>;
  InstrumentArrayType test_instrument_array_;
  const int n_instruments_{ test_helpers::RandomInt(2, 5) };
  void SetUp() override;
};

template <typename InputType>
void InstrumentationInstrumentArrayTest<InputType>::SetUp() {
  for (int i = 0; i < n_instruments_; ++i) {
    test_instrument_array_.AddInstrument(std::make_unique<InstrumentType>());
  }
}

using InstrumentTypes = ::testing::Types<moments::SphericalHarmonicI>;

TYPED_TEST_SUITE(InstrumentationInstrumentArrayTest, InstrumentTypes);

TYPED_TEST(InstrumentationInstrumentArrayTest, AddInstrument) {
  using InstrumentType = typename InstrumentationInstrumentArrayTest<TypeParam>::InstrumentType;
  using InstrumentArrayType = typename InstrumentationInstrumentArrayTest<TypeParam>::InstrumentArrayType;

  InstrumentArrayType test_instrument_array;
  const int n_instruments{ test_helpers::RandomInt(2, 5) };

  EXPECT_EQ(test_instrument_array.size(), 0);
  for (int i = 0; i < n_instruments; ++i) {
    EXPECT_NO_THROW({
      test_instrument_array.AddInstrument(std::make_unique<InstrumentType>());
    });
  }
  EXPECT_EQ(test_instrument_array.size(), n_instruments);
}

TYPED_TEST(InstrumentationInstrumentArrayTest, Iterators) {
  using InstrumentType = typename InstrumentationInstrumentArrayTest<TypeParam>::InstrumentType;
  for (auto& instrument : this->test_instrument_array_) {
    EXPECT_THAT(instrument.get(),
                WhenDynamicCastTo<InstrumentType*>(NotNull()));
  }
}

} // namespace
