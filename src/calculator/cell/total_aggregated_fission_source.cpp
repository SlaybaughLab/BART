#include "calculator/cell/total_aggregated_fission_source.hpp"

#include <deal.II/base/utilities.h>

#include "system/moments/spherical_harmonic_i.h"
#include "domain/domain_i.hpp"
#include "calculator/cell/integrated_fission_source_i.hpp"

namespace bart::calculator::cell {

template<int dim>
TotalAggregatedFissionSource<dim>::TotalAggregatedFissionSource(
    std::unique_ptr<IntegratedFissionSource> cell_fission_source_ptr,
    std::shared_ptr<Domain> domain_ptr)
    : cell_fission_source_ptr_(std::move(cell_fission_source_ptr)),
      domain_ptr_(domain_ptr) {
  this->set_description("Dealii aggregated fission source calculator", utility::DefaultImplementation(true));
  std::string call_location{"TotalAggregatedFissionSource constructor"};
  AssertPointerNotNull(cell_fission_source_ptr_.get(), "integrated fission source calculator", call_location);
  AssertPointerNotNull(domain_ptr_.get(), "domain pointer", call_location);
}


template<int dim>
auto TotalAggregatedFissionSource<dim>::AggregatedFissionSource(SystemMoments* system_moments_ptr) const -> double {
  auto cells = domain_ptr_->Cells();
  double fission_source{ 0 };
  for (auto& cell : cells) {
    fission_source += cell_fission_source_ptr_->CellValue(cell, system_moments_ptr);
  }
  return dealii::Utilities::MPI::sum(fission_source, MPI_COMM_WORLD);
}

template class TotalAggregatedFissionSource<1>;
template class TotalAggregatedFissionSource<2>;
template class TotalAggregatedFissionSource<3>;

} // namespace bart::calculator::cell