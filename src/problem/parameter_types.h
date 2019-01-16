#ifndef BART_SRC_PROBLEM_PARAMETER_TYPES_H_
#define BART_SRC_PROBLEM_PARAMETER_TYPES_H_

namespace bart {

namespace problem {

enum class AngularQuadType {
  kNone,
  kLevelSymmetricGaussChebyshev,
  kGaussLegendre,
};

enum class Boundary {
  kXMin = 0,
  kXMax = 1,
  kYMin = 2,
  kYMax = 3,
  kZMin = 4,
  kZMax = 5,
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

enum class InGroupSolverType {
  kNone,
  kSourceIteration,
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

enum class PreconditionerType {
  kNone,
  kAMG,
  kParaSails,
  kBlockJacobi,
  kJacobi,
  kSSOR,
};  

} // namespace problem

} // namespace bart

#endif // BART_SRC_PROBLEM_PARAMETER_TYPES_H_
