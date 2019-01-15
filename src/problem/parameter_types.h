#ifndef BART_SRC_PROBLEM_PARAMETER_TYPES_H_
#define BART_SRC_PROBLEM_PARAMETER_TYPES_H_

namespace bart {

namespace problem {

enum class AngularQuadType {
  kNone,
  kLevelSymmetricGaussChebyshev,
  kGaussLegendre,
};

enum class EigenSolverType {
  kNone,
  kPowerIteration,
};

enum class EquationType {
  kNone,
  kEvenParity,
  kSelfAdjointAngularFlux,
};

enum class LinearSolverType {
  kNone,
  kConjugateGradient,
  kGMRES,
  kBiCGSTAB,
  kDirect,
};

enum class MultiGroupSolverType {
  kNone,
  kGaussSeidel,
};

} // namespace problem

} // namespace bart

#endif // BART_SRC_PROBLEM_PARAMETER_TYPES_H_
