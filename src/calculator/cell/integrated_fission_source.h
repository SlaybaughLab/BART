#ifndef BART_SRC_CALCULATOR_CELL_INTEGRATED_FISSION_SOURCE_H_
#define BART_SRC_CALCULATOR_CELL_INTEGRATED_FISSION_SOURCE_H_

#include <memory>

#include "calculator/cell/integrated_fission_source_i.h"
#include "domain/domain_types.h"

namespace bart {

// Forward declarations of dependencies
namespace data {
struct CrossSections;
} // namespace data

namespace domain {
namespace finite_element {
template <int dim> class FiniteElementI;
} // namespace finite_element
} // namespace domain

namespace calculator {

namespace cell {

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
 *
 *
 * @tparam dim problem spatial dimension
 */
template <int dim>
class IntegratedFissionSource : public IntegratedFissionSourceI<dim> {
 public:
  /*! \brief Constructor, takes pointers to dependencies.
   *
   * This class takes shared ownership of a finite element object and a
   * cross-sections struct, both used to calculate the integrated cell fission
   * source.
   *
   * @param finite_element_ptr pointer to finite element object.
   * @param cross_sections_ptr pointer to cross-sections struct.
   */
  IntegratedFissionSource(
      std::shared_ptr<domain::finite_element::FiniteElementI<dim>> finite_element_ptr,
      std::shared_ptr<data::CrossSections> cross_sections_ptr);
  ~IntegratedFissionSource() = default;

  double CellValue(domain::CellPtr<dim> cell_ptr,
                   system::moments::SphericalHarmonicI* system_moments_ptr) const override;

 private:
  std::shared_ptr<domain::finite_element::FiniteElementI<dim>> finite_element_ptr_;
  std::shared_ptr<data::CrossSections> cross_sections_ptr_;
  const int cell_quadrature_points_;
};

} // namespace cell

} // namespace calculator

} // namespace bart

#endif // BART_SRC_CALCULATOR_CELL_INTEGRATED_FISSION_SOURCE_H_