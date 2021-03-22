#ifndef BART_SRC_CALCULATOR_CELL_INTEGRATED_FISSION_SOURCE_HPP_
#define BART_SRC_CALCULATOR_CELL_INTEGRATED_FISSION_SOURCE_HPP_

#include <memory>

#include "data/cross_sections/material_cross_sections.hpp"
#include "calculator/cell/integrated_fission_source_i.hpp"
#include "domain/domain_types.hpp"
#include "utility/has_dependencies.h"

namespace bart {

namespace domain::finite_element {
template <int dim> class FiniteElementI;
} // namespace domain::finite_element

namespace calculator::cell {

/*! \brief Calculates the integrated cell fission source.
 *
 * The total fission source is calculated using the following equation,
 * \f[
 *
 * f(\vec{r})= \sum_{g = 0}^G \nu_g\Sigma_{f,g}\phi_g(\vec{r})\;,
 *
 * \f]
 * where \f$G\f$ is the total number of groups. This class calculates the
 * integrated value using the cell quadrature provided by a
 * domain::finite_element::FiniteElementI object. It returns the value of
 *
 * \f[
 *
 * F = \int f(\vec{r})d\vec{r} = \sum_{q = 0}^Q f(\vec{r}_q)J(\vec{r}_q)
 * = \sum_{q = 0}^Q \sum_{g = 0}^G \nu_g \Sigma_{f,g} \phi_g(\vec{r}_q)  J(\vec{r}_q)
 *
 * \f]
 *
 * @tparam dim problem spatial dimension
 */
template <int dim>
class IntegratedFissionSource : public IntegratedFissionSourceI<dim>, public utility::HasDependencies {
 public:
  using CellPtr = domain::CellPtr<dim>;
  using MomentPtr = system::moments::SphericalHarmonicI*;
  using FiniteElement = domain::finite_element::FiniteElementI<dim>;
  /*! \brief Constructor, takes pointers to dependencies.
   *
   * This class takes shared ownership of a finite element object and a
   * cross-sections struct, both used to calculate the integrated cell fission
   * source.
   *
   * @param finite_element_ptr pointer to finite element object.
   * @param cross_sections_ptr pointer to cross-sections struct.
   */
  IntegratedFissionSource(std::shared_ptr<FiniteElement> finite_element_ptr,
                          std::shared_ptr<data::cross_sections::CrossSectionsI> cross_sections_ptr);
  ~IntegratedFissionSource() = default;

  data::cross_sections::CrossSectionsI* cross_sections_ptr() const { return cross_sections_ptr_.get(); };
  FiniteElement* finite_element_ptr() const { return finite_element_ptr_.get(); };

  [[nodiscard]] auto CellValue(CellPtr cell_ptr, MomentPtr system_moments_ptr) const -> double override;
 private:
  std::shared_ptr<FiniteElement> finite_element_ptr_;
  std::shared_ptr<data::cross_sections::CrossSectionsI> cross_sections_ptr_;
  const int cell_quadrature_points_;
};

} // namespace calculator::cell

} // namespace bart

#endif // BART_SRC_CALCULATOR_CELL_INTEGRATED_FISSION_SOURCE_HPP_