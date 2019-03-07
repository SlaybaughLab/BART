#ifndef BART_SRC_FORMULATION_EQUATION_TRANSPORT_I_H_
#define BART_SRC_FORMULATION_EQUATION_TRANSPORT_I_H_

#include <memory>

#include "formulation/types.h"
#include "domain/finite_element_i.h"


namespace bart {

namespace formulation {

namespace equation {

template <int dim>
class TransportI {
 public:
  using CellPtr = typename dealii::DoFHandler<dim>::active_cell_iterator;

  virtual ~TransportI() = default;

  /*! \brief Get the equation type (scalar/angular)
   *
   * \return formulation::EquationType that describes the type.
   */
  virtual EquationType equation_type() const = 0;

  /*! \brief Get the FEM discretization type (continuous/discontinuous)
   *
   * \return formulation::DiscretizationType that describes the type.
   */
  virtual DiscretizationType discretization_type() const = 0;

  virtual TransportI& ProvideFiniteElement(
      std::shared_ptr<domain::FiniteElementI<dim>>) = 0;

  /*! \brief Set the formulation to a specific cell.
   *
   * This is how cell-specific matrices will be generated.
   *
   * \param[in] to_set pointer to the cell to be set
   * \return this object
   */
  virtual TransportI& SetCell(CellPtr &to_set) = 0;
};

} // namespace formulation

} // namespace equation

} // namespace bart

#endif // BART_SRC_FORMULATION_EQUATION_TRANSPORT_H_