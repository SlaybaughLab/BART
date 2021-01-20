#ifndef BART_SRC_CALCULATOR_DRIFT_DIFFUSION_FACTORY_HPP_
#define BART_SRC_CALCULATOR_DRIFT_DIFFUSION_FACTORY_HPP_

#include "utility/factory/auto_registering_factory.hpp"

namespace bart::calculator::drift_diffusion {

template <int dim> class DriftDiffusionVectorCalculatorI;

enum class DriftDiffusionVectorCalculatorName {
  kDefaultImplementation = 0,
};

template <int dim, typename ...T>
class DriftDiffusionVectorCalculatorIFactory : public utility::factory::AutoRegisteringFactory<
    DriftDiffusionVectorCalculatorName,
    std::unique_ptr<DriftDiffusionVectorCalculatorI<dim>>(*)(T...)> {};

} // namespace bart::calculator::drift_diffusion

#endif //BART_SRC_CALCULATOR_DRIFT_DIFFUSION_FACTORY_HPP_
