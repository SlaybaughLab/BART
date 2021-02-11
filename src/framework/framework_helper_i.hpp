#ifndef BART_SRC_FRAMEWORK_FRAMEWORK_HELPER_I_HPP_
#define BART_SRC_FRAMEWORK_FRAMEWORK_HELPER_I_HPP_

#include "problem/parameters_i.h"
#include "framework/builder/framework_builder_i.hpp"
#include "material/material_protobuf.hpp"
#include "system/system_helper_i.hpp"
#include "system/moments/spherical_harmonic_i.h"

namespace bart::framework {

template <int dim>
class FrameworkHelperI {
 public:
  virtual ~FrameworkHelperI() = default;
  virtual auto ToFrameworkParameters(const problem::ParametersI& parameters) -> framework::FrameworkParameters = 0;
  virtual auto BuildFramework(builder::FrameworkBuilderI<dim>&,
                      framework::FrameworkParameters&) -> std::unique_ptr<framework::FrameworkI> = 0;
  virtual auto BuildFramework(builder::FrameworkBuilderI<dim>&,
                      framework::FrameworkParameters&,
                      system::moments::SphericalHarmonicI*) -> std::unique_ptr<framework::FrameworkI> = 0;
};


} // namespace bart::framework

#endif //BART_SRC_FRAMEWORK_FRAMEWORK_HELPER_I_HPP_
