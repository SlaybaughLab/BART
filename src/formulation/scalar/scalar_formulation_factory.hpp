#ifndef BART_SRC_FORMULATION_SCALAR_SCALAR_FORMULATION_FACTORY_HPP_
#define BART_SRC_FORMULATION_SCALAR_SCALAR_FORMULATION_FACTORY_HPP_

#include "utility/factory/auto_registering_factory.hpp"

namespace bart::formulation::scalar {

template <int dim> class DriftDiffusionI;

enum class DriftDiffusionFormulationName {
  kDefaultImplementation = 0,
};

template <int dim, typename ...T>
class DriftDiffusionIFactory : public utility::factory::AutoRegisteringFactory<
    DriftDiffusionFormulationName,
    std::unique_ptr<DriftDiffusionI<dim>>(*)(T...)> {};

} // namespace bart::formulation::scalar

#endif //BART_SRC_FORMULATION_SCALAR_SCALAR_FORMULATION_FACTORY_HPP_
