#ifndef BART_SRC_QUADRATURE_CALCULATORS_QUADRATURE_CALCULATORS_FACTORIES_HPP_
#define BART_SRC_QUADRATURE_CALCULATORS_QUADRATURE_CALCULATORS_FACTORIES_HPP_

#include "utility/factory/auto_registering_factory.hpp"

namespace bart::quadrature::calculators {

class AngularFluxIntegratorI;

enum class AngularFluxIntegratorName { kDefaultImplementation = 0 };

template <typename...T>
class AngularFluxIntegratorIFactory : public utility::factory::AutoRegisteringFactory<
    AngularFluxIntegratorName,
    std::unique_ptr<AngularFluxIntegratorI>(*)(T...)> {};

} // namespace bart::quadrature::calculators


#endif //BART_SRC_QUADRATURE_CALCULATORS_QUADRATURE_CALCULATORS_FACTORIES_HPP_
