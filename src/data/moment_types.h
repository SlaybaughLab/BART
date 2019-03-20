#ifndef BART_SRC_DATA_MOMENT_TYPES_H_
#define BART_SRC_DATA_MOMENT_TYPES_H_

#include <map>

#include <deal.II/lac/vector.h>

#include "utility/named_type.h"


namespace bart {

namespace data {

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

} // namespace data

} // namespace bart

#endif // BART_SRC_DATA_MOMENT_TYPES_H_