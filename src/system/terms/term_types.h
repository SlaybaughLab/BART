#ifndef BART_SRC_SYSTEM_TERMS_TERM_TYPES_H_
#define BART_SRC_SYSTEM_TERMS_TERM_TYPES_H_

#include "system/system_types.h"

namespace bart {

namespace system {

namespace terms {

//! Linear Terms that may vary iteration-to-iteration
enum class VariableLinearTerms {
  kOther = 0,            //!< Other source
  kScatteringSource = 1, //!< Scattering source
  kFissionSource = 2,    //!< Fission source
  kReflectiveBoundaryCondition = 3, //!< Reflective boundary conditions
};

//! Bilinear Terms that may vary iteration-to-iteration
enum class VariableBilinearTerms {
  kOther = 0,            //!< Other source
};

//! Standard pair types for terms
using MPILinearTermPair = std::pair<system::MPIVector, VariableLinearTerms>;
using MPIBilinearTermPair = std::pair<system::MPISparseMatrix, VariableBilinearTerms>;

} // namespace terms

} // namespace system

} // namespace bart

#endif // BART_SRC_SYSTEM_TERMS_TERM_TYPES_H_