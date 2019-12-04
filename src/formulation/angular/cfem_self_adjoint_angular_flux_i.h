#ifndef BART_SRC_FORMULATION_ANGULAR_CFEM_SELF_ADJOINT_ANGULAR_FLUX_I_H_
#define BART_SRC_FORMULATION_ANGULAR_CFEM_SELF_ADJOINT_ANGULAR_FLUX_I_H_

#include <deal.II/lac/full_matrix.h>
#include <deal.II/dofs/dof_accessor.h>

#include "formulation/formulation_types.h"
#include "quadrature/quadrature_point_i.h"
#include "system/system_types.h"

namespace bart {

namespace formulation {

namespace angular {

template <int dim>
class CFEMSelfAdjointAngularFluxI {
 public:
  struct InitializationToken{};

  virtual ~CFEMSelfAdjointAngularFluxI() = default;

  /*! \brief Integrates the bilinear streaming term and fills a given matrix.
   *
   * For a given cell in the triangulation, \f$K \in T_K\f$, with basis functions
   * \f$\varphi\f$, this function integrates bilinear SAAF streaming term for
   * one group using the cell quadrature adds them to the provided
   * local cell matrix, \f$\mathbf{A}\f$:
   * \f[
   * \mathbf{A}(i,j)_{K,g}' = \mathbf{A}(i,j)_{K,g} +
   * \int_{K}\left(\vec{\Omega}\cdot\nabla\varphi_i(\vec{r})\right)
   * \frac{1}{\sigma_{t,g}(\vec{r})}\left(\vec{\Omega}\cdot\nabla\varphi_j
   * (\vec{r})\right) dV
   * \f]
   *
   * where \f$\psi\f$ is the angular flux.
   *
   * @param to_fill cell matrix to fill.
   * @param init_token initialization token return by Initialize
   * @param cell_ptr pointer to the cell
   * @param quadrature_point quadrature angle to provide \f$\Omega\f$
   * @param group_number energy group to fill
   */
  virtual void FillCellStreamingTerm(
      FullMatrix& to_fill,
      const InitializationToken init_token,
      const CellPtr<dim>& cell_ptr,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const system::EnergyGroup group_number) const = 0;

  /*! \brief Initialize the formulation.
   * In general, this will pre-calculate matrix terms. The cell pointer is only
   * used to initialize the finite element object.
   * @param cell_ptr cell pointer for initialization.
   * @return Initialization token required for other calls
   */
  virtual InitializationToken Initialize(const formulation::CellPtr<dim>&) = 0;
};

} // namespace angular

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_ANGULAR_CFEM_SELF_ADJOINT_ANGULAR_FLUX_I_H_
