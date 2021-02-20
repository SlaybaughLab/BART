#ifndef BART_SRC_CALCULATOR_CELL_INTEGRATED_FISSION_SOURCE_I_HPP_
#define BART_SRC_CALCULATOR_CELL_INTEGRATED_FISSION_SOURCE_I_HPP_

#include "domain/domain_types.h"
#include "utility/has_description.h"

namespace bart {

namespace system::moments {
class SphericalHarmonicI;
} // namespace system::moments

namespace calculator::cell {

/*! \brief Interface for classes that calculate the cell-integrated fission source.
 *
 * This value is the scalar fission source for a single cell.
 *
 */
template <int dim>
class IntegratedFissionSourceI : public utility::HasDescription {
 public:
  using CellPtr = domain::CellPtr<dim>;
  using MomentPtr = system::moments::SphericalHarmonicI*;
  virtual ~IntegratedFissionSourceI() = default;

  /*! \brief Calculate and return the integrated cell fission source for a cell.
   *
   * @param cell_ptr pointer to the cell.
   * @param system_moments_ptr pointer to the system spherical harmonics.
   *
   * @return integrated cell fission source.
   */
  virtual auto CellValue(CellPtr cell_ptr, MomentPtr system_moments_ptr) const -> double = 0;
};

} // namespace calculator::cell

} // namespace bart

#endif // BART_SRC_CALCULATOR_CELL_INTEGRATED_FISSION_SOURCE_I_HPP_