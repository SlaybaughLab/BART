#ifndef BART_SRC_CALCULATOR_CELL_TOTAL_AGGREGATED_FISSION_SOURCE_I_H_
#define BART_SRC_CALCULATOR_CELL_TOTAL_AGGREGATED_FISSION_SOURCE_I_H_

namespace bart {

namespace system {

namespace moments {
class SphericalHarmonicI;
} // namespace moments

} // namespace system

namespace calculator {

namespace cell {

class TotalAggregatedFissionSourceI {
 public:
  virtual ~TotalAggregatedFissionSourceI() = default;
  virtual double AggregatedFissionSource(
      system::moments::SphericalHarmonicI* system_moments_ptr) const = 0;
};

} // namespace cell

} // namespace calculator

} // namespace bart

#endif // BART_SRC_CALCULATOR_CELL_TOTAL_AGGREGATED_FISSION_SOURCE_I_H_