#include "calculator/cell/total_aggregated_fission_source.h"

#include "system/moments/spherical_harmonic_i.h"
#include "domain/definition_i.h"
#include "calculator/cell/integrated_fission_source_i.h"

namespace bart {

namespace calculator {

namespace cell {

template<int dim>
double TotalAggregatedFissionSource<dim>::AggreatedFissionSource(
    system::moments::SphericalHarmonicI *system_moments_ptr) const {

  auto cells = domain_ptr_->Cells();
  double fission_source = 0;

  for (auto& cell : cells) {
    fission_source += cell_fission_source_ptr_->CellValue(cell, system_moments_ptr);
  }

  return fission_source;
}

template class TotalAggregatedFissionSource<1>;
template class TotalAggregatedFissionSource<2>;
template class TotalAggregatedFissionSource<3>;

} // namespace cell

} // namespace calculator

} // namespace bart