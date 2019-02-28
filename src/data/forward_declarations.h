#ifndef BART_SRC_DATA_FORWARD_DECLARATIONS_
#define BART_SRC_DATA_FORWARD_DECLARATIONS_

#include <map>
#include <memory>

#include <deal.II/lac/petsc_parallel_vector.h>

namespace bart {

namespace data {

struct SystemFluxes;

using Flux = dealii::PETScWrappers::MPI::Vector;
using MultiFluxPtrs = std::map<int, std::unique_ptr<Flux>>;

} // namespace data

} // namespace bart

#endif // BART_SRC_DATA_VECTOR_PARAMETERS_
