#include "matrix_parameters.h"

namespace bart {

namespace data {

std::shared_ptr<MPISparseMatrix> BuildMatrix(MatrixParameters &parameters) {
  auto matrix = std::make_shared<MPISparseMatrix>();
  matrix->reinit(parameters.rows, parameters.columns, parameters.dsp,
                 MPI_COMM_WORLD);
  return matrix;
}

} // namespace data

} // namespace bart
