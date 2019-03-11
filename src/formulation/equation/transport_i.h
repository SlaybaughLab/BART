#ifndef BART_SRC_FORMULATION_EQUATION_TRANSPORT_I_H_
#define BART_SRC_FORMULATION_EQUATION_TRANSPORT_I_H_

#include <memory>

#include "data/system_scalar_fluxes.h"
#include "data/forward_declarations.h"
#include "data/cross_sections.h"
#include "domain/finite_element_i.h"
#include "formulation/types.h"


namespace bart {

namespace formulation {

namespace equation {

template <int dim>
class TransportI {
 public:
  using CellPtr = typename dealii::DoFHandler<dim>::active_cell_iterator;
  using FaceNumber = int;

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

  /*! \brief Provide the finite element object.
   *
   * \param finite_element shared pointer to finite element object
   * \return this object
   */
  virtual TransportI& ProvideFiniteElement(
      std::shared_ptr<domain::FiniteElementI<dim>> finite_element) = 0;

  /*! \brief Provide the cross-section object.
   *
   * \param cross_sections shared pointer to cross-section object.
   * \return this object.
   */
  virtual TransportI& ProvideCrossSections(
      std::shared_ptr<data::CrossSections> cross_sections) = 0;

  virtual TransportI& ProvideScalarFluxes(
      std::shared_ptr<data::SystemScalarFluxes> scalar_fluxes) = 0;

  /*! \brief Set the formulation to a specific cell.
   *
   * This is how cell-specific matrices will be generated.
   *
   * \param[in] to_set pointer to the cell to be set
   */
  virtual void SetCell(const CellPtr &to_set) const = 0;

  virtual void SetFace(const CellPtr &to_set, FaceNumber face) const = 0;
};

} // namespace formulation

} // namespace equation

} // namespace bart

#endif // BART_SRC_FORMULATION_EQUATION_TRANSPORT_H_