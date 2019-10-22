#include "quadrature/factory/quadrature_factories.h"

#include "quadrature/ordinate.h"
#include "quadrature/quadrature_point.h"

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

template <int dim>
std::shared_ptr<QuadraturePointI<dim>> MakeQuadraturePointPtr(
    const QuadraturePointImpl impl) {
  std::shared_ptr<QuadraturePointI<dim>> quadrature_point_ptr = nullptr;

  if (impl == QuadraturePointImpl::kDefault) {
    // Default implementation
    quadrature_point_ptr = std::make_shared<QuadraturePoint<dim>>();
  }

  return quadrature_point_ptr;
}


template std::shared_ptr<OrdinateI<1>> MakeOrdinatePtr(const OrdinateType);
template std::shared_ptr<OrdinateI<2>> MakeOrdinatePtr(const OrdinateType);
template std::shared_ptr<OrdinateI<3>> MakeOrdinatePtr(const OrdinateType);

template std::shared_ptr<QuadraturePointI<1>> MakeQuadraturePointPtr(const QuadraturePointImpl);
template std::shared_ptr<QuadraturePointI<2>> MakeQuadraturePointPtr(const QuadraturePointImpl);
template std::shared_ptr<QuadraturePointI<3>> MakeQuadraturePointPtr(const QuadraturePointImpl);

} // namespace factory

} // namespace quadrature

} // namespace bart