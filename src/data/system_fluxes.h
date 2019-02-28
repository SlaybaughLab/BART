#ifndef BART_SRC_DATA_SYSTEM_FLUXES_H_
#define BART_SRC_DATA_SYSTEM_FLUXES_H_

#include <map>
#include <memory>

#include "data/forward_declarations.h"
#include "utility/uncopyable.h"

namespace bart {

namespace data {

struct SystemFluxes : private Uncopyable {
  MultiFluxPtrs current_iteration;
  MultiFluxPtrs previous_iteration;
};

} // namespace data

} // namespace bart

#endif // BART_SRC_DATA_SYSTEM_FLUXES_H_
