#ifndef BART_SRC_FORMULATION_TYPES_H_
#define BART_SRC_FORMULATION_TYPES_H_

#include <map>
#include <string>

namespace bart {

namespace formulation {

/*! Types of FEM Discretization available */
enum class DiscretizationType {
  kContinuousFEM = 0,    //!< Continuous FEM
  kDiscontinuousFEM = 1, //!< Discontinuous FEM
};

/*! Broad categories of formulations of the transport equation */
enum class EquationType {
  kScalar,  //!< Equations that solve for scalar flux directly
  kAngular, //!< Equations that solve for angular flux directly
};

/*! Specific scalar formulations of the transport equation */
enum class ScalarEquations {
  kDiffusion,
};

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_TYPES_H_