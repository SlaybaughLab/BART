#ifndef BART_SRC_PROBLEM_PARAMETER_TYPES_HPP_
#define BART_SRC_PROBLEM_PARAMETER_TYPES_HPP_

namespace bart {

namespace problem {

enum class AngularQuadType {
  kNone,
  kLevelSymmetricGaussian,
  kGaussLegendre
};

enum class DiscretizationType {
  kNone,
  kDiscontinuousFEM,
  kContinuousFEM,
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
  kDiffusion,
  kDriftDiffusion,
  kSelfAdjointAngularFlux,
  kTwoGridDiffusion
};

enum class CellFiniteElementType {
  kGaussian = 0,
};

enum class InGroupSolverType {
  kNone,
  kSourceIteration,
};

enum class LinearSolverType {
  kNone,
  kGMRES,
};

} // namespace problem

} // namespace bart

#endif // BART_SRC_PROBLEM_PARAMETER_TYPES_HPP_
