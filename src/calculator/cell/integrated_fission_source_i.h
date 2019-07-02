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

/*! \brief Interface for classes that calculate the cell norm of the fission
 *         source.
 *
 */

template <int dim>
class IntegratedFissionSourceI {
 public:
  virtual ~IntegratedFissionSourceI() = default;

  virtual double CellValue(domain::CellPtr<dim>,
                           system::moments::SphericalHarmonicI*) const = 0;
};

} // namespace cell

} // namespace calculator

} // namespace bart

#endif // BART_SRC_CALCULATOR_CELL_INTEGRATED_FISSION_SOURCE_I_H_