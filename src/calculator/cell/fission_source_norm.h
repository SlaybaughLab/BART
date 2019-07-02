#ifndef BART_SRC_CALCULATOR_CELL_INTEGRATED_FISSION_SOURCE_H_
#define BART_SRC_CALCULATOR_CELL_INTEGRATED_FISSION_SOURCE_H_

#include <memory>

#include "calculator/cell/fission_source_norm_i.h"
#include "domain/domain_types.h"

namespace bart {

// Forward declarations of dependencies
namespace data {
struct CrossSections;
} // namespace data

namespace domain {
template <int dim> class FiniteElementI;
} // namespace domain

namespace calculator {

namespace cell {

/*! \brief Calculates the cell norm of the fission source.
 *
 * The total fission source is calculated using the following equation,
 * \f[
 *
 * f(\vec{r})= \sum_{g = 0}^G \nu_g\Sigma_{f,g}\phi_g(\vec{r})\;,
 *
 * \f]
 * where \f$G\f$ is the total number of groups.
 *
 *
 *
 * @tparam dim problem spatial dimension
 */
template <int dim>
class IntegratedFissionSource : public IntegratedFissionSourceI<dim> {
 public:
  IntegratedFissionSource(
      std::shared_ptr<domain::FiniteElementI<dim>> finite_element_ptr,
      std::shared_ptr<data::CrossSections> cross_sections_ptr);
  ~IntegratedFissionSource() = default;

  double GetCellNorm(domain::CellPtr<dim> cell_ptr,
                     system::moments::SphericalHarmonicI* system_moments_ptr) const override;

 private:
  std::shared_ptr<domain::FiniteElementI<dim>> finite_element_ptr_;
  std::shared_ptr<data::CrossSections> cross_sections_ptr_;
  const int cell_quadrature_points_;
};

} // namespace cell

} // namespace calculator

} // namespace bart

#endif // BART_SRC_CALCULATOR_CELL_INTEGRATED_FISSION_SOURCE_H_