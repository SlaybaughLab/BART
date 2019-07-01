#ifndef BART_SRC_CALCULATOR_CELL_FISSION_SOURCE_NORM_H_
#define BART_SRC_CALCULATOR_CELL_FISSION_SOURCE_NORM_H_

#include "calculator/cell/fission_source_norm_i.h"

namespace bart {

namespace calculator {

namespace cell {

/*! \brief Calculates the cell norm of the fission source.
 *
 * @tparam dim problem spatial dimension
 */
template <int dim>
class FissionSourceNorm : public FissionSourceNormI {};

} // namespace cell

} // namespace calculator

} // namespace bart

#endif // BART_SRC_CALCULATOR_CELL_FISSION_SOURCE_NORM_H_