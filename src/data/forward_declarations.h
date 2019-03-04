#ifndef BART_SRC_DATA_FORWARD_DECLARATIONS_
#define BART_SRC_DATA_FORWARD_DECLARATIONS_

#include <map>
#include <memory>

#include <deal.II/lac/petsc_parallel_vector.h>

namespace bart {

namespace data {

// Variable types
using GroupNumber = int;   //!< Problem group number.
using OrdinateIndex = int; //!< Angular quadrature ordinate index
using HarmonicL = int;     //!< Spherical harmonic \f$\ell\f$ value for moments
using HarmonicM = int;     //!< Spherical harmonic \f$m\f$ value for moments

// System flux types
using Flux = dealii::PETScWrappers::MPI::Vector;
using ScalarFluxPtrs = std::map<GroupNumber, std::unique_ptr<Flux>>;

// Flux data structures
struct CrossSections;
struct SystemScalarFluxes;

// Flux builders



} // namespace data

} // namespace bart

#endif // BART_SRC_DATA_VECTOR_PARAMETERS_
