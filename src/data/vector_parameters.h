#ifndef BART_SRC_DATA_VECTOR_PARAMETERS_
#define BART_SRC_DATA_VECTOR_PARAMETERS_

#include <map>

#include <deal.II/lac/petsc_parallel_vector.h>

namespace bart {

namespace data {

typedef dealii::PETScWrappers::MPI::Vector Flux;
typedef std::map<int, Flux> GroupFluxes;

} // namespace data

} // namespace bart

#endif // BART_SRC_DATA_VECTOR_PARAMETERS_
