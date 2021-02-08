#ifndef BART_SRC_CALCULATOR_CELL_TESTS_TOTAL_AGGREGATED_FISSION_SOURCE_MOCK_HPP_
#define BART_SRC_CALCULATOR_CELL_TESTS_TOTAL_AGGREGATED_FISSION_SOURCE_MOCK_HPP_

#include "calculator/cell/total_aggregated_fission_source_i.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace system {
namespace moments {
class SphericalHarmonicI;
} // namespace moments
} // namespace system

namespace calculator {

namespace cell {

class TotalAggregatedFissionSourceMock : public TotalAggregatedFissionSourceI {
 public:
  MOCK_CONST_METHOD1(AggregatedFissionSource, double(
      system::moments::SphericalHarmonicI* system_moments_ptr));
};

} // namespace cell

} // namespace calculator

} // namespace bart

#endif // BART_SRC_CALCULATOR_CELL_TESTS_TOTAL_AGGREGATED_FISSION_SOURCE_MOCK_HPP_