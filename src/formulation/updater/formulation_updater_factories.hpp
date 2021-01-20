#ifndef BART_SRC_FORMULATION_UPDATER_FORMULATION_UPDATER_FACTORIES_HPP_
#define BART_SRC_FORMULATION_UPDATER_FORMULATION_UPDATER_FACTORIES_HPP_

#include "utility/factory/auto_registering_factory.hpp"

namespace bart::formulation::updater {

template <int dim> class DriftDiffusionUpdater;

enum class DriftDiffusionUpdaterName {
  kDefaultImplementation = 0,
};

template <int dim, typename ...T>
class DriftDiffusionUpdaterFactory : public utility::factory::AutoRegisteringFactory<
    DriftDiffusionUpdaterName,
    std::unique_ptr<DriftDiffusionUpdater<dim>>(*)(T...)> {};

} // namespace bart::formulation::updater

#endif //BART_SRC_FORMULATION_UPDATER_FORMULATION_UPDATER_FACTORIES_HPP_
