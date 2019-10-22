#include "quadrature/factory/quadrature_factories.h"

#include "quadrature/ordinate.h"

namespace bart {

namespace quadrature {

namespace factory {

template <int dim>
std::shared_ptr<OrdinateI<dim>> MakeOrdinatePtr(const OrdinateType type) {

  std::shared_ptr<OrdinateI<dim>> ordinate_ptr = nullptr;

  if (type == OrdinateType::kDefault) {
    // Default implementation
    ordinate_ptr = std::make_shared<Ordinate<dim>>();
  }

  return ordinate_ptr;
}

template std::shared_ptr<OrdinateI<1>> MakeOrdinatePtr<1>(const OrdinateType);
template std::shared_ptr<OrdinateI<2>> MakeOrdinatePtr<2>(const OrdinateType);
template std::shared_ptr<OrdinateI<3>> MakeOrdinatePtr<3>(const OrdinateType);

} // namespace factory

} // namespace quadrature

} // namespace bart