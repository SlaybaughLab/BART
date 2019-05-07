#ifndef BART_SRC_DATA_SYSTEM_RHS_LHS_TYPES_H_
#define BART_SRC_DATA_SYSTEM_RHS_LHS_TYPES_H_

#include <array>
#include <utility>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>

namespace bart {

namespace data {

namespace system {

//! Group number for rhs and lhs
using GroupNumber = int;
//! Angle index for rhs and lhs
using AngleIndex = int;
//! Index used to store and access rhs vectors and lhs matrices
using Index = std::pair<GroupNumber, AngleIndex>;

//! Sparse MPI matrix used for left-hand-side matrices
using MPISparseMatrix = dealii::PETScWrappers::MPI::SparseMatrix;
//! Sparse MPI vector used for ride-hand-side vectors
using MPIVector = dealii::PETScWrappers::MPI::Vector;

using HarmonicL = int;
//!< Spherical harmonic \f$\ell\f$ value for moments

using HarmonicM = int;
//!< Spherical harmonic \f$m\f$ value for moments

/*! \typedef HarmonicIndex
 * \brief Index of a spherical harmonic in the form \f$[g, \ell, m]\f$.
 *
 */
using MomentIndex = std::array<int, 3>;

using MomentVector = dealii::Vector<double>; //!< Vector for storing moments

using MomentsMap = std::map<MomentIndex, MomentVector>;

} // namespace system

} // namespace data

} // namespace bart

#endif // BART_SRC_DATA_SYSTEM_RHS_LHS_TYPES_H_