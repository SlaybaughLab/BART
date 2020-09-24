#include "instrumentation/factory/instrument_factories.h"

#include <sstream>
#include <deal.II/lac/vector.h>

#include "instrumentation/instrument.h"
#include "instrumentation/outstream/to_ostream.h"
#include "instrumentation/converter/multi_converter.h"
#include "instrumentation/converter/calculator/vector_subtractor.h"
#include "instrumentation/converter/fourier/fourier_transform.h"
#include "instrumentation/converter/dealii_to_complex_vector.h"
#include "instrumentation/converter/to_string/int_vector_complex_pair_to_string.h"
#include "instrumentation/outstream/to_ostream.h"

#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

namespace instrumentation = bart::instrumentation;
namespace factory = bart::instrumentation::factory;
namespace test_helpers = bart::test_helpers;
using ::testing::WhenDynamicCastTo, ::testing::NotNull;

TEST(InstrumentationFactoryInstrumentFactoriesTest,
     MakeErrorFourierTransformInstrument) {
  /* The Error fourier transform instrument should string together three
   * converters and an outstream. We will verify all the correct components have
   * been built and in the right order. */
  // Check that the pointer returned is of type Instrument that reads a dealiivector
  using DealiiVector = dealii::Vector<double>;

  const int vector_size{test_helpers::RandomInt(10, 20)};
  DealiiVector error_vector(vector_size);
  for (int i = 0; i < vector_size; ++i) {
    error_vector[i] = test_helpers::RandomDouble(-100, 100);
  }

  std::unique_ptr<std::ostream> string_stream = std::make_unique<std::ostringstream>();
  auto instrument_ptr =
      factory::MakeErrorFourierTransformInstrument(error_vector,
                                                   std::move(string_stream));
  ASSERT_NE(instrument_ptr, nullptr);
}

} // namespace
