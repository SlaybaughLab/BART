#ifndef BART_SRC_DATA_SYSTEM_SCALAR_FLUXES_H_
#define BART_SRC_DATA_SYSTEM_SCALAR_FLUXES_H_

#include <map>
#include <memory>

#include "data/forward_declarations.h"
#include "utility/uncopyable.h"

namespace bart {

namespace data {

struct SystemScalarFluxes : private utility::Uncopyable {
  data::ScalarFluxPtrs current_iteration;
  data::ScalarFluxPtrs previous_iteration;
};

std::unique_ptr<data::ScalarFluxPtrs> BuildSystemScalarFluxes(
    int total_groups,
    dealii::IndexSet locally_owned_dofs);

} // namespace data

} // namespace bart

#endif // BART_SRC_DATA_SYSTEM_SCALAR_FLUXES_H_
