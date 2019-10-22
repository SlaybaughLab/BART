#ifndef BART_SRC_QUADRATURE_FACTORY_QUADRATURE_FACTORIES_H_
#define BART_SRC_QUADRATURE_FACTORY_QUADRATURE_FACTORIES_H_

#include <memory>

#include "quadrature/ordinate_i.h"
#include "quadrature/quadrature_point_i.h"
#include "quadrature/quadrature_types.h"

namespace bart {

namespace quadrature {

namespace factory {

template <int dim>
std::shared_ptr<OrdinateI<dim>> MakeOrdinatePtr(
    const OrdinateType type = OrdinateType::kDefault);

template <int dim>
std::shared_ptr<QuadraturePointI<dim>> MakeQuadraturePointPtr(
    const QuadraturePointImpl impl = QuadraturePointImpl::kDefault);

} // namespace factory

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_FACTORY_QUADRATURE_FACTORIES_H_
