#ifndef BART_SRC_FORMULATION_SCALAR_DIFFUSION_I_HPP_
#define BART_SRC_FORMULATION_SCALAR_DIFFUSION_I_HPP_

#include <deal.II/lac/full_matrix.h>
#include <deal.II/dofs/dof_accessor.h>

#include "system/moments/spherical_harmonic_types.h"
#include "utility/has_description.h"

//! Scalar (non-angular) formulations of the transport equation.
namespace bart::formulation::scalar {

/*! \brief Interface for classes that provide the diffusion formulation of the transport equation.
 *
 * The diffusion formulation is a common non-angular version of the transport equation that models neutron propegation
 * diffusively, much like the heat equation. Therefore, this is naturally a second-order formulation of the transport
 * equation. The diffusion equation is derived by integrating the first-order transport equation over angle and
 * using Fick's law as a closure for the current. The multigroup diffusion equation for a multiplying medium
 * in _k_-eigenvalue form is,
 *
 * \f[
 * -D_g(\vec{r})\vec{\nabla} \cdot \vec{\nabla} \phi_g(\vec{r}) + (\Sigma_{t,g} - \Sigma_{s}^{g\to g})\phi_g(\vec{r})
 * = \sum_{g' \neq g}\Sigma_{g}^{g' \to g}\phi_{g'}(\vec{r}) + \frac{\chi_g}{k}\sum_{g'}\nu_{g'}\Sigma_{f,g'}\phi_{g'}(\vec{r})
 * \f]
 *
 * With weak formulation,
 *
 * \f[
 * \bigg(\nabla \cdot v(\vec{r}), D_g\nabla\cdot\phi_{g}(\vec{r})\bigg)_{K} + \bigg(v(\vec{r}), (\Sigma_{t,g} - \Sigma_{s}^{g\to g})\phi_g(\vec{r})\bigg)_{K}
 * + \bigg(v(\vec{r}), \frac{1}{2}\phi_g(\vec{r})\bigg)_{\partial K, \text{vacuum}} =
 * \bigg(v(\vec{r}),\sum_{g' \neq g}\Sigma_{g}^{g' \to g}\phi_{g'}(\vec{r})\bigg)_{K} +
 * \bigg(v(\vec{r}), \frac{\chi_g}{k}\sum_{g'}\nu_{g'}\Sigma_{f,g'}\phi_{g'}(\vec{r})\bigg)_{K}\;.
 * \f]
 *
 * Each of these terms are stamped individually on the system matrix, \f$\mathbf{A}\f$ or vector \f$\vec{b}\f$
 * by the member functions of classes that derive from this interface. The integration over the element \f$K\f$ is done
 * using the cell quadrature and the right-hand-side scalar fluxes are treated as sources from the previous iteration.
 *
 * For further information about this derivation see any neutron transport text such as
 * <a href="https://www.ans.org/store/item-350016/">Lewis and Miller</a>.
 * @tparam dim spatial dimension
 */
template <int dim>
class DiffusionI : public utility::HasDescription {
 public:
  //! Types of boundaries for the diffusion equation
  enum class BoundaryType {
    kVacuum,
    kReflective
  };

  //! Pointer to a cell iterator returned by a dof object.
  using CellPtr = typename dealii::DoFHandler<dim>::active_cell_iterator;
  using Matrix = dealii::FullMatrix<double>;
  using Vector = dealii::Vector<double>;
  using GroupNumber = int;
  using FaceNumber = int;

  virtual ~DiffusionI() = default;

  /*! \brief Precalculate many of the shape functions.
   *
   * The bilinear terms require the shape-funtion or gradient of the shape function squared. This function precalculates
   * those, reducing the number of times the underlying finite element object needs to be called. The shape functions
   * for each cell are identical, and the jacobian is used to translate from the base cell.
   *
   * @param cell_ptr an arbitrary cell to use, often just the beginning of active cells.
   */
  virtual auto Precalculate(const CellPtr& cell_ptr) -> void = 0;

  /*! \brief Integrates the bilinear streaming term over a cell and fills a given matrix.
   *
   * For a given cell in the triangulation, \f$K \in T_K\f$, with basis functions
   * \f$\varphi\f$, this function integrates the bilinear SAAF streaming term for
   * one group using the cell quadrature and adds them to the provided
   * local cell matrix, \f$\mathbf{A}\f$:
   *
   * \f[
   * \mathbf{A}(i,j)_{K,g}' = \mathbf{A}(i,j)_{K,g} + \int_{K}D_g\nabla\varphi_i(\vec{r})\nabla\varphi_j(\vec{r})dV
   * \f]
   *
   * @param to_fill
   */
  virtual auto FillCellStreamingTerm(Matrix& to_fill, const CellPtr&, GroupNumber) const -> void = 0;
  /*! \brief Integrates the bilinear collision term over a cell and fills a given matrix.
   *
   * For a given cell in the triangulation, \f$K \in T_K\f$, with basis functions
   * \f$\varphi\f$, this function integrates the bilinear SAAF streaming term for
   * one group using the cell quadrature and adds them to the provided
   * local cell matrix, \f$\mathbf{A}\f$:
   *
   * \f[
   * \mathbf{A}(i,j)_{K,g}' = \mathbf{A}(i,j)_{K,g} + \int_{K}\varphi_i(\vec{r})(\Sigma_{t,g} - \Sigma_{s, g \to g})\varphi_j(\vec{r})dV
   * \f]
   *
   * @param to_fill
   */
  virtual auto FillCellCollisionTerm(Matrix& to_fill, const CellPtr&, GroupNumber) const -> void = 0;
/*! \brief Integrates the bilinear boundary term over a cell and fills a given matrix.
   *
   * For a given cell and face in the triangulation, \f$\partial K \in \partial T_K\f$, with basis functions
   * \f$\varphi\f$, this function integrates the bilinear diffusion boundary term for
   * one group using the cell quadrature and adds them to the provided
   * local cell matrix, \f$\mathbf{A}\f$. For reflective boundary conditions, nothing is added, for vacuum:
   *
   * \f[
   * \mathbf{A}(i,j)_{K,g}' = \mathbf{A}(i,j)_{K,g} + \frac{1}{2}\int_{\partial K}\varphi_i(\vec{r})\varphi_j(\vec{r})dV
   * \f]
   *
   * @param to_fill
   */
  virtual auto FillBoundaryTerm(Matrix& to_fill, const CellPtr&, FaceNumber, BoundaryType) const -> void = 0;

  /*! \brief Integrates the fixed source term and fills a given vector.
   *
   * For a given cell in the triangulation, \f$K \in T_K\f$, with basis functions
   * \f$\varphi\f$, this function integrates the diffusion fixed-source term and
   * adds it to the cell right-hand side vector.
   * \f[
   * \vec{b}(i)_{K,g}' = \vec{b}(i)_{K,g} + \int_{K}q_{g}(\vec{r})\varphi_i(\vec{r})dV
   * \f]
   *
   * where \f$\phi\f$ is the scalar flux. Adds the result per cell DOFF to the
   * input-output vector cell_rhs.
   *
   *
   * @param to_fill cell vector to fill
   */
  virtual auto FillCellFixedSource(Vector& to_fill, const CellPtr&, GroupNumber) const -> void = 0;

  /*! \brief Integrates the fission source term and fills a given cell vector.
   *
   * For a given cell in the triangulation, \f$K \in T_K\f$, with basis functions
   * \f$\varphi\f$, this function integrates the diffusion fission-source term and
   * adds it to the cell right-hand side vector.
   * \f[
   * \vec{b}(i)_{K,g}' = \vec{b}(i)_{K,g} + \frac{\chi_g}{k}\sum_{g'}\int_{K}\varphi_i(\vec{r})\Sigma_{f,g'}\nu_{g'}\phi_{g'}(\vec{r})dV
   * \f]
   *
   * @param to_fill cell vector to fill
   * @param k_eigenvalue value of the _k_-eigenvalue
   * @param in_group_moment scalar flux moment for the current group
   * @param group_moments scalar flux moment for all other groups
   */
  virtual auto FillCellFissionSource(Vector& to_fill, const CellPtr&, GroupNumber, double k_eigenvalue,
                                     const system::moments::MomentVector& in_group_moment,
                                     const system::moments::MomentsMap& group_moments) const -> void = 0;

  /*! \brief Integrates the scattering source term and fills a given cell vector.
   *
   * For a given cell in the triangulation, \f$K \in T_K\f$, with basis functions
   * \f$\varphi\f$, this function integrates the diffusion scattering-source term and
   * adds it to the cell right-hand side vector.
   * \f[
   * \vec{b}(i)_{K,g}' = \vec{b}(i)_{K,g} + \sum_{g' \neq g}\int_{K}\varphi_i(\vec{r})\Sigma_{s,g' \to g}\phi_{g'}(\vec{r})dV
   * \f]
   *
   * @param to_fill cell vector to fill
   * @param k_eigenvalue value of the _k_-eigenvalue
   * @param group_moments scalar flux moment for all groups
   */
  virtual auto FillCellScatteringSource(Vector& to_fill, const CellPtr&, GroupNumber,
                                        const system::moments::MomentsMap& group_moments) const -> void = 0;

  /*! \brief Returns a bool indicating if Precalculate been called. */
  virtual auto is_initialized() const -> bool = 0;
};

} // namespace bart::formulation::scalar

#endif //BART_SRC_FORMULATION_SCALAR_DIFFUSION_I_HPP_