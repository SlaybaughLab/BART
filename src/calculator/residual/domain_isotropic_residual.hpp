#ifndef BART_SRC_CALCULATOR_RESIDUAL_DOMAIN_ISOTROPIC_RESIDUAL_HPP_
#define BART_SRC_CALCULATOR_RESIDUAL_DOMAIN_ISOTROPIC_RESIDUAL_TEST_HPP_

#include "calculator/residual/domain_isotropic_residual_i.hpp"
#include "domain/domain_i.hpp"
#include "calculator/residual/cell_isotropic_residual_i.hpp"
#include "utility/has_dependencies.h"

namespace bart::calculator::residual {

/*! \brief Default implementation of domain isotropic scattering residual calculator. */
template <int dim>
 class DomainIsotropicResidual : public DomainIsotropicResidualI, public utility::HasDependencies {
 public:
   using Domain = typename domain::DomainI<dim>;
   using CellIsotropicResidualCalculator = typename calculator::residual::CellIsotropicResidualI<dim>;
   using Vector = dealii::Vector<double>;

   /*! \brief Constructor.
    * Takes ownership of a CellIsotropicResidualCalculator and shared ownership of a domain.
    */
   DomainIsotropicResidual(std::unique_ptr<CellIsotropicResidualCalculator>, std::shared_ptr<Domain>);

   auto CalculateDomainResidual(FluxMoments *current_flux_moments,
                                FluxMoments *previous_flux_moments) -> Vector override;

   //! Access cell isotropic residual calculator dependency.
   auto cell_isotropic_residual_calculator_ptr() { return cell_isotropic_residual_calculator_ptr_.get(); }
   //! Access domain dependency.
   auto domain_ptr() { return domain_ptr_.get(); }
 private:
   std::unique_ptr<CellIsotropicResidualCalculator> cell_isotropic_residual_calculator_ptr_;
   std::shared_ptr<Domain> domain_ptr_;
};

} // namespace bart::calculator::residual

#endif //BART_SRC_CALCULATOR_RESIDUAL_DOMAIN_ISOTROPIC_RESIDUAL_HPP_
