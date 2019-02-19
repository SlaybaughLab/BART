#ifndef BART_SRC_DATA_MATRIX_PARAMETERS_
#define BART_SRC_DATA_MATRIX_PARAMETERS_

#include <deal.II/base/index_set.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

/*! This struct holds the parameters necessary for creation of a parallel
 * PETSc matrix. */

namespace bart {

namespace data {

struct MatrixParameters {
  dealii::IndexSet rows;
  dealii::IndexSet columns;
  dealii::DynamicSparsityPattern dsp;  
};

} // namespace data

} // namespace bart

#endif // BART_SRC_DATA_MATRIX_PARAMETERS_
