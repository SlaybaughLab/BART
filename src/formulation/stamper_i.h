#ifndef BART_SRC_FORMULATION_STAMPER_I_H_
#define BART_SRC_FORMULATION_STAMPER_I_H_

#include <functional>

#include "domain/domain_types.h"
#include "formulation/formulation_types.h"
#include "system/system_types.h"

namespace bart {

namespace formulation {

/*! \brief Interface for a system matrix stamper.
 * See the default implementation for detailed description for the purpose of
 * this class.
 *
 * \tparam dim spatial dimension of the cells in the mesh.
 *
 * \author J.S. Rehak
 */
template <int dim>
class StamperI {
 public:
  virtual ~StamperI() = default;
  virtual void StampMatrix(
      system::MPISparseMatrix& to_stamp,
      std::function<void(formulation::FullMatrix&,
                         const domain::CellPtr<dim>&)> stamping_function) = 0;
};

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_STAMPER_I_H_