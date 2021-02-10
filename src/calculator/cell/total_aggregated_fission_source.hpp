#ifndef BART_SRC_CALCULATOR_CELL_TOTAL_AGGREGATED_FISSION_SOURCE_HPP_
#define BART_SRC_CALCULATOR_CELL_TOTAL_AGGREGATED_FISSION_SOURCE_HPP_

#include <memory>

#include "calculator/cell/total_aggregated_fission_source_i.hpp"
#include "calculator/cell/integrated_fission_source_i.hpp"
#include "domain/domain_types.h"

namespace bart {

namespace domain {
template <int dim> class DefinitionI;
} // namespace domain

namespace calculator::cell {

template <int dim>
class TotalAggregatedFissionSource : public TotalAggregatedFissionSourceI {
 public:
  using SystemMoments = system::moments::SphericalHarmonicI;

  TotalAggregatedFissionSource(std::unique_ptr<IntegratedFissionSourceI<dim>> cell_fission_source_ptr,
                               std::shared_ptr<domain::DefinitionI<dim>> domain_ptr)
      : cell_fission_source_ptr_(std::move(cell_fission_source_ptr)),
        domain_ptr_(domain_ptr) {};
  virtual ~TotalAggregatedFissionSource() = default;

  [[nodiscard]] auto AggregatedFissionSource(SystemMoments*) const -> double override;

  IntegratedFissionSourceI<dim>* cell_fission_source_ptr() const { return cell_fission_source_ptr_.get(); }

 private:
  std::unique_ptr<IntegratedFissionSourceI<dim>> cell_fission_source_ptr_;
  std::shared_ptr<domain::DefinitionI<dim>> domain_ptr_;
};

} // namespace calculator::cell

} // namespace bart

#endif // BART_SRC_CALCULATOR_CELL_TOTAL_AGGREGATED_FISSION_SOURCE_HPP_