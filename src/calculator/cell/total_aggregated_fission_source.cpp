#include "calculator/cell/total_aggregated_fission_source.hpp"

#include <deal.II/base/utilities.h>

#include "system/moments/spherical_harmonic_i.h"
#include "domain/definition_i.h"
#include "calculator/cell/integrated_fission_source_i.hpp"

namespace bart::calculator::cell {

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