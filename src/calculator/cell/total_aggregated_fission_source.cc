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
  return 0;
}

template class TotalAggregatedFissionSource<1>;
template class TotalAggregatedFissionSource<2>;
template class TotalAggregatedFissionSource<3>;

} // namespace cell

} // namespace calculator

} // namespace bart