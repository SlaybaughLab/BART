#ifndef BART_SRC_CALCULATOR_CELL_TESTS_TOTAL_AGGREGATED_FISSION_SOURCE_MOCK_HPP_
#define BART_SRC_CALCULATOR_CELL_TESTS_TOTAL_AGGREGATED_FISSION_SOURCE_MOCK_HPP_

#include "calculator/cell/total_aggregated_fission_source_i.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace system::moments {
class SphericalHarmonicI;
} // namespace system::moments

namespace calculator::cell {

class TotalAggregatedFissionSourceMock : public TotalAggregatedFissionSourceI {
 public:
  MOCK_METHOD(double, AggregatedFissionSource, (system::moments::SphericalHarmonicI*), (const, override));
};

} // namespace calculator::cell

} // namespace bart

#endif // BART_SRC_CALCULATOR_CELL_TESTS_TOTAL_AGGREGATED_FISSION_SOURCE_MOCK_HPP_