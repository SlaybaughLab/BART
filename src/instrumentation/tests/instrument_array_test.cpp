#include "instrumentation/instrument_array.hpp"

#include <deal.II/lac/vector.h>

#include "convergence/status.hpp"
#include "instrumentation/tests/instrument_mock.h"
#include "system/moments/spherical_harmonic_i.h"
#include "system/moments/tests/spherical_harmonic_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"
#include "utility/colors.h"

namespace  {

namespace convergence = bart::convergence;
namespace instrumentation = bart::instrumentation;
namespace moments = bart::system::moments;
namespace utility = bart::utility;
namespace test_helpers = bart::test_helpers;

using ::testing::WhenDynamicCastTo, ::testing::NotNull, ::testing::Ref;

template <typename InputType>
class InstrumentationInstrumentArrayTest : public ::testing::Test {
 public:
  using InstrumentType = instrumentation::InstrumentMock<InputType>;
  using InstrumentArrayType = instrumentation::InstrumentArray<InputType>;
  InstrumentArrayType test_instrument_array_;

  // Test parameters
  const int n_instruments_{ test_helpers::RandomInt(2, 5) };
  std::unique_ptr<InputType> test_input_;

  // Helper functions
  std::unique_ptr<InputType> GetTestInput() const;

  void SetUp() override;
};

template <typename InputType>
void InstrumentationInstrumentArrayTest<InputType>::SetUp() {
  for (int i = 0; i < n_instruments_; ++i) {
    test_instrument_array_.AddInstrument(std::make_unique<InstrumentType>());
  }
  test_input_ = GetTestInput();
}

template <>
auto InstrumentationInstrumentArrayTest<moments::SphericalHarmonicI>::GetTestInput()  const
-> std::unique_ptr<moments::SphericalHarmonicI> {
  return std::make_unique<moments::SphericalHarmonicMock>();
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

TYPED_TEST(InstrumentationInstrumentArrayTest, Read) {
  using InstrumentType = typename InstrumentationInstrumentArrayTest<TypeParam>::InstrumentType;
  auto& test_input = *this->test_input_;

  for (auto& instrument : this->test_instrument_array_) {
    auto dynamic_ptr = dynamic_cast<InstrumentType*>(instrument.get());
    EXPECT_CALL(*dynamic_ptr, Read(Ref(test_input)));
  }

  this->test_instrument_array_.Read(test_input);
}

} // namespace
