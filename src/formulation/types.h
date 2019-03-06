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

/*! Maps discetization types to strings */
std::map<DiscretizationType, std::string> discretization_type_name_map {
    {DiscretizationType::kContinuousFEM, "continuous FEM"},
    {DiscretizationType::kDiscontinuousFEM, "discontinuous FEM"}
};

/*! Broad categories of formulations of the transport equation */
enum class EquationType {
  kScalar,  //!< Equations that solve for scalar flux directly
  kAngular, //!< Equations that solve for angular flux directly
};

/*! Maps equation types to strings */
std::map<EquationType, std::string> equation_type_name_map {
    {EquationType::kScalar, "scalar"},
    {EquationType::kAngular, "angular"}
};

/*! Specific scalar formulations of the transport equation */
enum class ScalarEquations {
  kDiffusion,
};

std::map<ScalarEquations, std::string> scalar_equations_name_map {
    {ScalarEquations::kDiffusion, "diffusion equation"},
};

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_TYPES_H_