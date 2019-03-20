#ifndef BART_SRC_DATA_FORWARD_DECLARATIONS_
#define BART_SRC_DATA_FORWARD_DECLARATIONS_

#include <map>
#include <memory>

#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>

namespace bart {

namespace data {

// Variable types
using GroupNumber = int;   //!< Problem group number.
using OrdinateIndex = int; //!< Angular quadrature ordinate index



// System flux types
using FluxVector = dealii::PETScWrappers::MPI::Vector;
using FluxVectorPtr = std::unique_ptr<data::FluxVector>;
using RightHandSideVector = dealii::PETScWrappers::MPI::Vector;

using ScalarRightHandSidePtrs =
    std::map<GroupNumber, std::unique_ptr<RightHandSideVector>>;
using ScalarFluxPtrs = std::map<GroupNumber, std::unique_ptr<FluxVector>>;

// System Matrices
using SystemMatrix = dealii::PETScWrappers::MPI::SparseMatrix;
using ScalarSystemMatrixPtrs = std::map<GroupNumber, std::unique_ptr<SystemMatrix>>;

// Flux data structures
struct CrossSections;
struct SystemScalarFluxes;

// Flux builders



} // namespace data

} // namespace bart

#endif // BART_SRC_DATA_VECTOR_PARAMETERS_
