#include "data/system_scalar_fluxes.h"

namespace bart {

namespace data {

std::shared_ptr<data::ScalarFluxPtrs> BuildSystemScalarFluxes(
    int total_groups,
    dealii::IndexSet locally_owned_dofs) {

  auto flux_ptrs = std::make_shared<data::ScalarFluxPtrs>();

  for (int group = 0; group < total_groups; ++group) {

    auto flux = std::make_unique<data::Flux>();
    flux->reinit(locally_owned_dofs, MPI_COMM_WORLD);

    flux_ptrs->insert(std::make_pair(group, std::move(flux)));
  }

  return std::move(flux_ptrs);
}

} // namespace data

} // namespace bart


