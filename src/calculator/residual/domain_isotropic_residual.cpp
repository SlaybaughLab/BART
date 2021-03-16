#include "calculator/residual/domain_isotropic_residual.hpp"

namespace bart::calculator::residual {

template <int dim>
DomainIsotropicResidual<dim>::DomainIsotropicResidual(
    std::unique_ptr<CellIsotropicResidualCalculator> cell_isotropic_residual_calculator_ptr,
    std::shared_ptr<Domain> domain_ptr)
    : cell_isotropic_residual_calculator_ptr_(std::move(cell_isotropic_residual_calculator_ptr)),
      domain_ptr_(std::move(domain_ptr)) {
  std::string function_name{ "DomainIsotropicResidual constructor"};
  this->AssertPointerNotNull(domain_ptr_.get(), "domain", function_name);
  this->AssertPointerNotNull(cell_isotropic_residual_calculator_ptr_.get(), "cell isotropic residual calculator",
                             function_name);
}

template class DomainIsotropicResidual<1>;
template class DomainIsotropicResidual<2>;
template class DomainIsotropicResidual<3>;

} // namespace bart::calculator::residual
