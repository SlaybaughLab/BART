#ifndef BART_SRC_DATA_VECTOR_PARAMETERS_
#define BART_SRC_DATA_VECTOR_PARAMETERS_

#include <deal.II/lac/petsc_parallel_vector.h>

namespace bart {

namespace data {

typedef dealii::PETScWrappers::MPI::Vector Flux;

} // namespace data

} // namespace bart

#endif // BART_SRC_DATA_VECTOR_PARAMETERS_
