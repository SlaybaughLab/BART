#ifndef BART_SRC_PROBLEM_PARAMETER_TYPES_H_
#define BART_SRC_PROBLEM_PARAMETER_TYPES_H_

namespace bart {

namespace problem {

enum class AngularQuadType {
  kLevelSymmetricGaussChebyshev,
  kGaussLegendre,
  kNone,
};

enum class EigenSolverType {
  kPowerIteration,
};

enum class EquationType {
  kEvenParity,
  kSelfAdjointAngularFlux,
  kNone
};

enum class LinearSolverType {
  kConjugateGradient,
  kGMRES,
  kBiCGSTAB,
  kDirect,
};

enum class MultiGroupSolverType {
  kGaussSeidel,
};

} // namespace problem

} // namespace bart

#endif // BART_SRC_PROBLEM_PARAMETER_TYPES_H_
