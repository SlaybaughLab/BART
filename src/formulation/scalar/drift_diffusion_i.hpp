#ifndef BART_SRC_FORMULATION_SCALAR_DRIFT_DIFFUSION_I_HPP_
#define BART_SRC_FORMULATION_SCALAR_DRIFT_DIFFUSION_I_HPP_

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <quadrature/quadrature_types.h>

#include "domain/domain_types.hpp"
#include "formulation/formulation_types.hpp"
#include "system/system_types.h"

namespace bart::formulation::scalar {

/*! \brief Interface for classes that provide the drift diffusion formulation terms.
 *
 * Classes the derive from this provide two terms from the drift-diffusion formulation used by the non-linear diffusion
 * acceleration method (<a href="https://doi.org/10.13182/NSE11-81">Park, et. al. 2012</a>). The drift-diffusion
 * formulation for the _k_-eigenvalue problem is:
 *
 * \f[
 * -D_g(\vec{r})\vec{\nabla} \cdot \vec{\nabla} \phi_g(\vec{r}) - \vec{D}_g(\vec{r})\vec{\nabla}\phi_g(\vec{r}) +
 * (\Sigma_{t,g} - \Sigma_{s}^{g\to g})\phi_g(\vec{r})
 * = \sum_{g' \neq g}\Sigma_{g}^{g' \to g}\phi_{g'}(\vec{r}) + \frac{\chi_g}{k}\sum_{g'}\nu_{g'}\Sigma_{f,g'}\phi_{g'}(\vec{r})
 * \f]
 *
 * Where \f$\vec{D}_g\f$ is the drift-diffusion term. There are multiple ways to describe this term, and this class is
 * not expected to calculate it directly, but using a class of type calculator::drift_diffusion::DriftDiffusionVectorCalculatorI.
 * This calculation requires the current, so the function to fill that term requires the current at each quadrature
 * point.
 *
 * @tparam dim
 */
template <int dim>
class DriftDiffusionI {
 public:
  using BoundaryType = formulation::BoundaryType;
  using CellPtr = typename domain::CellPtr<dim>;
  using EnergyGroup = system::EnergyGroup;
  using Matrix = typename dealii::FullMatrix<double>;
  using Vector = typename dealii::Vector<double>;
  using VectorMap = std::map<quadrature::QuadraturePointIndex, std::shared_ptr<Vector>>;
  virtual ~DriftDiffusionI() = default;

  /*! \brief Integrates the bilinear boundary term over a cell and fills a given matrix.
   *
   * For a given cell and face in the triangulation, \f$\partial K \in \partial T_K\f$, with basis functions
   * \f$\varphi\f$, this function integrates the bilinear drift-diffusion boundary term for
   * one group using the cell quadrature and adds them to the provided
   * local cell matrix
   *
   * @param to_fill cell matrix to fill
   * @param group_angular_flux angular flux for all groups
   */
  virtual auto FillCellBoundaryTerm(Matrix& to_fill, const CellPtr&, domain::FaceIndex, BoundaryType,
                                    const VectorMap& group_angular_flux) const -> void = 0;

  /*! \brief Integrates the bilinear drift-diffusion term over a cell and fills a given matrix.
   *
   * For a given cell in the triangulation, \f$K \in T_K\f$, with basis functions
   * \f$\varphi\f$, this function integrates the bilinear SAAF streaming term for
   * one group using the cell quadrature and adds them to the provided
   * local cell matrix, \f$\mathbf{A}\f$.
   *
   * \f[
   * \mathbf{A}(i,j)_{K,g}' = \mathbf{A}(i,j)_{K,g} + \int_{K}\vec{D}_g\varphi_j(\vec{r})\nabla\varphi_i(\vec{r})dV
   * \f]
   *
   * @param to_fill
   */
  virtual auto FillCellDriftDiffusionTerm(Matrix& to_fill, const CellPtr&, system::EnergyGroup,
                                          const Vector& group_scalar_flux,
                                          const std::array<Vector, dim>& current) const -> void = 0;
};

} // namespace bart::formulation::scalar

#endif //BART_SRC_FORMULATION_SCALAR_DRIFT_DIFFUSION_I_HPP_
