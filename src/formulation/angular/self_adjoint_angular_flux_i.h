#ifndef BART_SRC_FORMULATION_ANGULAR_SELF_ADJOINT_ANGULAR_FLUX_I_H_
#define BART_SRC_FORMULATION_ANGULAR_SELF_ADJOINT_ANGULAR_FLUX_I_H_

#include <deal.II/lac/full_matrix.h>
#include <deal.II/dofs/dof_accessor.h>

#include "domain/domain_types.hpp"
#include "formulation/formulation_types.hpp"
#include "quadrature/quadrature_point_i.hpp"
#include "system/system_types.h"
#include "system/moments/spherical_harmonic_types.h"

namespace bart {

namespace formulation {

namespace angular {

template <int dim>
class SelfAdjointAngularFluxI {
 public:
  virtual ~SelfAdjointAngularFluxI() = default;

  /*!
  * \brief Integrates the bilinear boundary terms and fills a given matrix.
  *
  * For a given boundary in the triangulation, \f$\partial K \in \partial T_K\f$,
  * with basis functions \f$\varphi\f$, this function integrates the bilinear
  * boundary term for one group using the quadrature specified in the
  * problem definition and adds them to the local cell matrix. For the case
  * \f$(\vec{n} \cdot \vec{\Omega}) > 0\f$, the boundary condition is
  * \f$\Psi_b(\vec{r},\vec{\Omega}) = \Psi(\vec{r},\vec{\Omega})\f$ and is
  * partially handled by integrating the following bilinear term (there may
  * also be a linear term):
  * \f[
  * \mathbf{A}(i,j)_{K,g}' = \mathbf{A}(i,j)_{K,g} +
  * \int_{\partial K}
  * (\hat{n}\cdot\vec{\Omega})\varphi_i(\vec{r})
  * \varphi_j(\vec{r})
  * dS
  * \f]
  *
  * \return No values returned, modifies input parameter \f$\mathbf{A}\to \mathbf{A}'\f$.
  !*/
  virtual void FillBoundaryBilinearTerm(
      FullMatrix& to_fill,
      const domain::CellPtr<dim>& cell_ptr,
      const domain::FaceIndex face_number,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const system::EnergyGroup group_number) = 0;

  /*! \brief Fills the linear boundary term for reflective boundary conditions.
   *
   */
   virtual auto FillReflectiveBoundaryLinearTerm(
      Vector& to_fill,
      const domain::CellPtr<dim>& cell_ptr,
      const domain::FaceIndex face_number,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const dealii::Vector<double>& incoming_flux) -> double = 0;

  /*!
   * \brief Integrates the bilinear collision term and fills a given matrix.
   *
   * For a given cell in the triangulation, \f$K \in T_K\f$, with basis functions
   * \f$\varphi\f$, this function integrates the bilinear SAAF collision terms
   * for one group using the cell quadrature and adds them to the local cell
   * matrix.
   * \f[
   * \mathbf{A}(i,j)_{K,g}' = \mathbf{A}(i,j)_{K,g} +
   * \int_{K}\sigma_{t,g}(\vec{r})\varphi_i(\vec{r})
   * \varphi_j(\vec{r}) dV
   * \f]
   *
   * where \f$\psi\f$ is the angular flux.
   *
   * @param to_fill cell matrix to fill.
   * @param cell_ptr pointer to the cell
   * @param group_number energy group to fill
   * \return No values returned, modifies input parameter \f$\mathbf{A}\to \mathbf{A}'\f$.
   */
  virtual void FillCellCollisionTerm(
      FullMatrix& to_fill,
      const domain::CellPtr<dim>& cell_ptr,
      const system::EnergyGroup group_number) = 0;

  /*! \brief Integrates the linear fission terms and fills a given vector.
   *
   * For a given cell in the triangulation, \f$K \in T_K\f$, with basis functions
   * \f$\varphi\f$, this function integrates the two fission source SAAF terms
   * for one group using the quadrature specified in the problem definition and
   * adds them to the cell right-hand side vector.
   *
   * \f[
   * \vec{b}(i)_{K,g}' = \vec{b}(i)_{K,g} +
   * \int_{K}\frac{q_g(\vec{r})}{4\pi}\varphi_i(\vec{r})dV +
   * \int_{K}\left(\vec{\Omega}\cdot\nabla\varphi_i(\vec{r})\right)
   * \frac{q_g(\vec{r})}{4\pi\sigma_{t,g}(\vec{r})}dV
   * \f]
   *
   * where \f$q_g(\vec{r})\f$ is the fission source:
   * \f[
   * q_g(\vec{r}) = \frac{\nu_g \cdot \sigma_{f,g}(\vec{r})_g}
   * {4\pi k_{\text{eff}}}\phi_g(\vec{r})
   * \f]
   *
   *
   * @param to_fill cell vector to fill
   * @param cell_ptr pointer to the cell
   * @param quadrature_point quadrature point
   * @param group_number energy group number
   * @param k_eff k effective
   * @param in_group_moment in-group flux moments
   * @param group_moments full set of group moments
   */
  virtual auto FillCellFissionSourceTerm(
      Vector& to_fill,
      const domain::CellPtr<dim>& cell_ptr,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const system::EnergyGroup group_number,
      const double k_eff,
      const system::moments::MomentVector& in_group_moment,
      const system::moments::MomentsMap& group_moments) -> double = 0;

  /*!
 * \brief Integrates the linear fixed-source terms and fills a given vector.
 *
 * For a given cell in the triangulation, \f$K \in T_K\f$, with basis functions
 * \f$\varphi\f$, this function integrates the two linear fixed source SAAF terms
 * for one group using the quadrature specified in the problem definition and
 * adds them to the cell right-hand side vector.
 *
 * \f[
 * \vec{b}(i)_{K,g}' = \vec{b}(i)_{K,g} +
 * \int_{K}\frac{q_g(\vec{r})}{4\pi}\varphi_i(\vec{r})dV +
 * \int_{K}\left(\vec{\Omega}\cdot\nabla\varphi_i(\vec{r})\right)
 * \frac{q_g(\vec{r})}{4\pi\sigma_{t,g}(\vec{r})}dV
 * \f]
 *
 * where \f$q_g(\vec{r})\f$ is a given fixed source.
 */
  virtual void FillCellFixedSourceTerm(
      Vector& to_fill,
      const domain::CellPtr<dim>& cell_ptr,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const system::EnergyGroup group_number) = 0;

  /*! \brief Integrates the linear scattering source and fills a given vector.
   *
   * For a given cell in the triangulation, \f$K \in T_K\f$, with basis functions
   * \f$\varphi\f$, this function integrates the two linear scattering SAAF terms
   * for one group using the quadrature specified in the problem definition and
   * adds them to the cell right-hand side vector.
   * \f[
   * \vec{b}(i)_{K,g}' = \vec{b}(i)_{K,g} + \sum_{g'}\left[
   * \int_{K}\frac{\sigma_{s,g'\to g}(\vec{r})}{4\pi}\phi_{g'}(\vec{r}) \varphi_i(\vec{r}) dV +
   * \int_{K}\left(\vec{\Omega}\cdot \nabla \varphi_i(\vec{r})\right)
   * \frac{\sigma_{s,g'\to g}(\vec{r})}{4\pi\sigma_{t,g}(\vec{r}) }\phi_{g'}(\vec{r}) dV
   * \right]
   * \f]
   *
   * where \f$\phi\f$ is the scalar flux. Adds the result per cell DOFF to the
   * input-output vector cell_rhs.
   *
   *
   * @param to_fill cell vector to fill
   * @param cell_ptr pointer to the cell
   * @param quadrature_point quadrature point to provide \f$\Omega\f$
   * @param in_group_moment in-group scalar flux moment
   * @param group_moments out-group scalar flux moments
   */
  virtual auto FillCellScatteringSourceTerm(
      Vector& to_fill,
      const domain::CellPtr<dim>& cell_ptr,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const system::EnergyGroup group_number,
      const system::moments::MomentVector& in_group_moment,
      const system::moments::MomentsMap& group_moments) -> double = 0;

  /*! \brief Integrates the bilinear streaming term and fills a given matrix.
   *
   * For a given cell in the triangulation, \f$K \in T_K\f$, with basis functions
   * \f$\varphi\f$, this function integrates the bilinear SAAF streaming term for
   * one group using the cell quadrature and adds them to the provided
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
   * @param cell_ptr pointer to the cell
   * @param quadrature_point quadrature angle to provide \f$\Omega\f$
   * @param group_number energy group to fill
   * \return No values returned, modifies input parameter \f$\mathbf{A}\to \mathbf{A}'\f$.
   */
  virtual void FillCellStreamingTerm(
      FullMatrix& to_fill,
      const domain::CellPtr<dim>& cell_ptr,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const system::EnergyGroup group_number) = 0;


  /*! \brief Initialize the formulation.
   * In general, this will pre-calculate matrix terms. The cell pointer is only
   * used to initialize the finite element object.
   * @param cell_ptr cell pointer for initialization.
   */
  virtual void Initialize(const domain::CellPtr<dim>&) = 0;
};

} // namespace angular

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_ANGULAR_SELF_ADJOINT_ANGULAR_FLUX_I_H_
