#ifndef BART_SRC_CALCULATOR_CELL_TOTAL_AGGREGATED_FISSION_SOURCE_I_HPP_
#define BART_SRC_CALCULATOR_CELL_TOTAL_AGGREGATED_FISSION_SOURCE_I_HPP_

namespace bart {

namespace system::moments {
class SphericalHarmonicI;
} // namespace system::moments

namespace calculator::cell {

class TotalAggregatedFissionSourceI {
 public:
  virtual ~TotalAggregatedFissionSourceI() = default;
  virtual auto AggregatedFissionSource(system::moments::SphericalHarmonicI* system_moments_ptr) const -> double = 0;
};

} // namespace calculator::cell

} // namespace bart

#endif // BART_SRC_CALCULATOR_CELL_TOTAL_AGGREGATED_FISSION_SOURCE_I_HPP_