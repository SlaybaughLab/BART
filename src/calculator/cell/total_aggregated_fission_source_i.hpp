#ifndef BART_SRC_CALCULATOR_CELL_TOTAL_AGGREGATED_FISSION_SOURCE_I_HPP_
#define BART_SRC_CALCULATOR_CELL_TOTAL_AGGREGATED_FISSION_SOURCE_I_HPP_

namespace bart {

namespace system::moments {
class SphericalHarmonicI;
} // namespace system::moments

//! Calculators for various cell quantities.
namespace calculator::cell {
/*! \brief Interface for classes that calculate the total aggregated fission source for a domain.
 *
 * This is distinct from the integrated fission source over a cell, in that it returns the total for the whole domain.
 *
 */
class TotalAggregatedFissionSourceI {
 public:
  virtual ~TotalAggregatedFissionSourceI() = default;
  /*! \brief Calculate the aggregated fission source for the domain.
   *
   * @param system_moments_ptr system moments to be used to calculate the fission source.
   * @return double value of the domain fission source.
   */
  virtual auto AggregatedFissionSource(system::moments::SphericalHarmonicI* system_moments_ptr) const -> double = 0;
};

} // namespace calculator::cell

} // namespace bart

#endif // BART_SRC_CALCULATOR_CELL_TOTAL_AGGREGATED_FISSION_SOURCE_I_HPP_