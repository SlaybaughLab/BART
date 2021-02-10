#ifndef BART_SRC_CALCULATOR_CELL_TOTAL_AGGREGATED_FISSION_SOURCE_HPP_
#define BART_SRC_CALCULATOR_CELL_TOTAL_AGGREGATED_FISSION_SOURCE_HPP_

#include <memory>

#include "calculator/cell/total_aggregated_fission_source_i.hpp"
#include "calculator/cell/integrated_fission_source_i.hpp"
#include "domain/domain_types.h"
#include "utility/has_dependencies.h"
#include "utility/has_description.h"

namespace bart {

namespace domain {
template <int dim> class DefinitionI;
} // namespace domain

namespace calculator::cell {

template <int dim>
class TotalAggregatedFissionSource : public TotalAggregatedFissionSourceI, public utility::HasDescription,
                                     public utility::HasDependencies {
 public:
  using Domain = domain::DefinitionI<dim>;
  using IntegratedFissionSource = IntegratedFissionSourceI<dim>;
  using SystemMoments = system::moments::SphericalHarmonicI;

  TotalAggregatedFissionSource(std::unique_ptr<IntegratedFissionSource>, std::shared_ptr<Domain>);
  virtual ~TotalAggregatedFissionSource() = default;

  [[nodiscard]] auto AggregatedFissionSource(SystemMoments*) const -> double override;

  IntegratedFissionSource* cell_fission_source_ptr() const { return cell_fission_source_ptr_.get(); }
  Domain* domain_ptr() const { return domain_ptr_.get(); }
 private:
  std::unique_ptr<IntegratedFissionSource> cell_fission_source_ptr_;
  std::shared_ptr<Domain> domain_ptr_;
};

} // namespace calculator::cell

} // namespace bart

#endif // BART_SRC_CALCULATOR_CELL_TOTAL_AGGREGATED_FISSION_SOURCE_HPP_