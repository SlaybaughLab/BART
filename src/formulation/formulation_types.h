#ifndef BART_SRC_FORMULATION_FORMULATION_TYPES_H_
#define BART_SRC_FORMULATION_FORMULATION_TYPES_H_

namespace bart {

namespace formulation {

using FullMatrix = dealii::FullMatrix<double>;
using Vector = dealii::Vector<double>;

enum class BoundaryType {
  kVacuum = 0,
  kReflective = 1
};

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_FORMULATION_TYPES_H_
