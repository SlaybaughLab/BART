#ifndef BART_SRC_CALCULATOR_CELL_TOTAL_AGGREGATED_FISSION_SOURCE_H_
#define BART_SRC_CALCULATOR_CELL_TOTAL_AGGREGATED_FISSION_SOURCE_H_

#include "domain/domain_types.h"

#include <memory>

#include "calculator/cell/total_aggregated_fission_source_i.h"

namespace bart {

namespace calculator {

namespace cell {
template <int dim>
class IntegratedFissionSourceI;

template <int dim>
class TotalAggregatedFissionSource : public TotalAggregatedFissionSourceI<dim> {
 public:
  TotalAggregatedFissionSource(
      std::unique_ptr<IntegratedFissionSourceI<dim>> cell_fission_source_ptr)
      : cell_fission_source_ptr_(std::move(cell_fission_source_ptr)) {};
  virtual ~TotalAggregatedFissionSource() = default;

 private:
  std::unique_ptr<IntegratedFissionSourceI<dim>> cell_fission_source_ptr_;

};

} // namespace cell

} // namespace calculator

} // namespace bart

#endif // BART_SRC_CALCULATOR_CELL_TOTAL_AGGREGATED_FISSION_SOURCE_H_