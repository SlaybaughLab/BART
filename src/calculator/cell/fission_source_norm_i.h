#ifndef BART_SRC_CALCULATOR_CELL_FISSION_SOURCE_NORM_I_H_
#define BART_SRC_CALCULATOR_CELL_FISSION_SOURCE_NORM_I_H_

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
 */

template <int dim>
class FissionSourceNormI {
 public:
  virtual ~FissionSourceNormI() = default;

  virtual double GetCellNorm(domain::CellPtr<dim> cell_ptr,
                             system::moments::SphericalHarmonicI* system_moments) const = 0;
};

} // namespace cell

} // namespace calculator

} // namespace bart

#endif // BART_SRC_CALCULATOR_CELL_FISSION_SOURCE_NORM_I_H_