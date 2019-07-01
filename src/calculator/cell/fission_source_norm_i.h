#ifndef BART_SRC_CALCULATOR_CELL_FISSION_SOURCE_NORM_I_H_
#define BART_SRC_CALCULATOR_CELL_FISSION_SOURCE_NORM_I_H_

namespace bart {

namespace calculator {

namespace cell {

/*! \brief Interface for classes that calculate the cell norm of the fission
 *         source.
 */

class FissionSourceNormI {
 public:
  virtual ~FissionSourceNormI() = default;
};

} // namespace cell

} // namespace calculator

} // namespace bart

#endif // BART_SRC_CALCULATOR_CELL_FISSION_SOURCE_NORM_I_H_