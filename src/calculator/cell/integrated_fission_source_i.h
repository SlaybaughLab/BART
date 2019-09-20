#ifndef BART_SRC_CALCULATOR_CELL_INTEGRATED_FISSION_SOURCE_I_H_
#define BART_SRC_CALCULATOR_CELL_INTEGRATED_FISSION_SOURCE_I_H_

#include "domain/domain_types.h"

namespace bart {

namespace system {
namespace moments {
class SphericalHarmonicI;
} // namespace moments
} // namespace system

namespace calculator {

namespace cell {

/*! \brief Interface for classes that calculate the integrated cell fission
 *         source.
 *
 */

template <int dim>
class IntegratedFissionSourceI {
 public:
  virtual ~IntegratedFissionSourceI() = default;

  /*! \brief Calculate and return the integrated cell fission source for a cell.
   *
   * @param cell_ptr pointer to the cell.
   * @param system_moments_ptr pointer to the system spherical harmonics.
   *
   * @return integrated cell fission source.
   */
  virtual double CellValue(
      domain::CellPtr<dim> cell_ptr,
      system::moments::SphericalHarmonicI* system_moments_ptr) const = 0;
};

} // namespace cell

} // namespace calculator

} // namespace bart

#endif // BART_SRC_CALCULATOR_CELL_INTEGRATED_FISSION_SOURCE_I_H_