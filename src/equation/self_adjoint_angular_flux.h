#ifndef BART_SRC_EQUATION_SELF_ADJOINT_ANGULAR_FLUX_H_
#define BART_SRC_EQUATION_SELF_ADJOINT_ANGULAR_FLUX_H_

#include <algorithm>

#include "equation_base.h"

/*!
 * \brief Provides the weak formulation for the Self Adjoint Angular Flux 
 *         (SAAF) equation.
 *
 * SAAF is a second-order formulation of the transport equation, given as:
 * \f[
 * -\vec{\Omega} \cdot \nabla \frac{1}{\sigma_t}\vec{\Omega}\cdot\nabla\psi
 * + \sigma_t\psi = \frac{1}{4\pi}\left[\sigma_s\phi + Q - \vec{\Omega}\cdot
 * \nabla\frac{\sigma_s\phi + Q}{\sigma_t}\right]
 * \f]
 * Where for eigenvalue problems:
 *
 * \f[
 * Q = 
 * \frac{\nu\sigma_f}{4\pi k_{\text{eff}}}\phi
 * \f]
 *
 * 
 * and with boundary conditions:
 * \f[
 * \psi_b = \cases{
 * F(\vec{r}, \vec{\Omega}) & for $\hat{n} \cdot \vec{\Omega} < 0$ \\
 * \psi(\vec{r}, \vec{\Omega}) & for $\hat{n} \cdot \vec{\Omega} > 0$
 * }
 * \f]
 *
 * This class provides a Continuous Galerkin formulation that includes either
 * vacuum boundary conditions or reflective boundary conditions.
 * 
 * \author Joshua Rehak
 */

template <int dim>
class SelfAdjointAngularFlux : public EquationBase<dim> {
 public:
  /*! Class constructor.
   * \param equation_name name of the equation
   * \param prm ParameterHandler object containing problem definition
   * \param data_ptr pointer to FundamentalData holding problem data
   */
  SelfAdjointAngularFlux(const std::string equation_name,
                         const dealii::ParameterHandler &prm,
                         std::shared_ptr<FundamentalData<dim>> &data_ptr);

  /*! Default class destructor */
  ~SelfAdjointAngularFlux() = default;

  void AssembleLinearForms (const int &g) override;
  
  /*!
   * \brief Integrates the bilinear boundary terms in the SAAF equation and adds
   *        the values to the matrix cell_matrix.
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
   * \param cell the cell \f$K\f$.
   * \param fn face index for the current face in the cell
   * \param cell_matrix the local matrix to be modified, \f$\mathbf{A}\f$
   * \param g energy group
   * \param dir direction \f$\vec{\Omega}\f$.
   * \return No values returned, modifies input parameter \f$\mathbf{A}\to \mathbf{A}'\f$.   
   !*/
  void IntegrateBoundaryBilinearForm (
      typename dealii::DoFHandler<dim>::active_cell_iterator &,
      const int &,
      dealii::FullMatrix<double> &cell_matrix,
      const int &,
      const int &dir) override;
  /*!
   * \brief Integrates the linear bilinear boundary term in the SAAF equation and adds
   *        the values to the matrix cell_matrix.
   *           
   * For a given boundary in the triangulation, \f$\partial K \in \partial T_K\f$,
   * with basis functions \f$\varphi\f$, this function integrates the linear
   * boundary term for one group using the quadrature specified in the
   * problem definition and adds them to the local cell right-hand side. For the case
   * \f$(\vec{n} \cdot \vec{\Omega}) < 0\f$, the boundary condition is
   * \f$\Psi_b(\vec{r},\vec{\Omega}) = \Psi_{\text{inc}}(\vec{r},\vec{\Omega})\f$ and is
   * partially handled by integrating the following bilinear term (there is also
   * a bilinear term):
   * \f[
   * \vec{b}(i)_{K,g}' = \vec{b}(i)_{K,g} +
   * \int_{\partial K}
   * (\hat{n}\cdot\vec{\Omega})\varphi_i(\vec{r})
   * \Psi_{\text{inc},g}(\vec{r})
   * dS
   * \f]
   *
   * For vacuum boundary conditions, \f$\Psi_{\text{inc},g} = 0\f$. For reflective
   * boundary conditions, it is equal to the angular flux on the boundary in
   * the previous iteration.
   * 
   * \param cell the cell \f$K\f$.
   * \param fn face index for the current face in the cell
   * \param cell_rhs the local cell right-hand side to be modified, \f$\vec{b}\f$
   * \param g energy group
   * \param dir direction \f$\vec{\Omega}\f$.
   * \return No values returned, modifies input parameter \f$\vec{b}\to \vec{b}'\f$.   
   !*/
  void IntegrateBoundaryLinearForm (
      typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
      const int &fn,/*face number*/
      dealii::Vector<double> &cell_rhs,
      const int &g,
      const int &dir) override;  
  
  /*!
   * \brief Integrates the bilinear terms in the SAAF equation and adds
   *        the values to the matrix cell_matrix.
   *           
   * For a given cell in the triangulation, \f$K \in T_K\f$, with basis functions
   * \f$\varphi\f$, this function integrates two bilinear SAAF terms for one group
   * using the quadrature specified in the problem definition and adds them to
   * the local cell matrix.
   * \f[
   * \mathbf{A}(i,j)_{K,g}' = \mathbf{A}(i,j)_{K,g} + 
   * \int_{K}\left(\vec{\Omega}\cdot\nabla\varphi_i(\vec{r})\right)
   * \frac{1}{\sigma_{t,g}(\vec{r})}\left(\vec{\Omega}\cdot\nabla\varphi_j
   * (\vec{r})\right) dV + \int_{K}\sigma_{t,g}(\vec{r})\varphi_i(\vec{r})
   * \varphi_j(\vec{r}) dV
   * \f]
   *
   * where \f$\psi\f$ is the angular flux. 
   *
   * \param cell the cell \f$K\f$.
   * \param cell_matrix the local matrix to be modified, \f$\mathbf{A}\f$
   * \param g energy group
   * \param dir direction \f$\vec{\Omega}\f$.
   * \return No values returned, modifies input parameter \f$\mathbf{A}\to \mathbf{A}'\f$.
   */
  void IntegrateCellBilinearForm (
      typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
      dealii::FullMatrix<double> &cell_matrix,
      const int &g,
      const int &dir) override;

  /*!
   * \brief Integrates the linear fixed terms in the SAAF equation and adds the
   *        values to the vector cell_rhs.
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
   * where \f$q_g(\vec{r})\f$ is a given fixed source or a fission source:
   *
   * \f[
   * q_g(\vec{r}) = \frac{\nu_g \cdot \sigma_{f,g}(\vec{r})_g}
   * {4\pi k_{\text{eff}}}\phi_g(\vec{r})
   * \f]
   */   
  void IntegrateCellFixedLinearForm (
      typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
      dealii::Vector<double> &cell_rhs,
      const int &g,
      const int &dir) override;
  
  /*!
   * \brief Integrates the linear scattering terms in the SAAF equation and adds
   *        the values to the vector cell_rhs.
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
   * \param cell the cell \f$K\f$.
   * \param cell_rhs the local right-hand side to be modified, \f$\vec{b}\f$
   * \param g energy group
   * \param dir direction \f$\vec{\Omega}\f$.
   * \return No values returned, modifies input parameter \f$\vec{b}\to \vec{b}'\f$.
   */
  void IntegrateScatteringLinearForm (
      typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
      dealii::Vector<double> &cell_rhs,
      const int &g,
      const int &dir) override;
  
  void PreassembleCellMatrices () override;



  void IntegrateInterfaceBilinearForm (
      typename dealii::DoFHandler<dim>::active_cell_iterator &,
      typename dealii::DoFHandler<dim>::cell_iterator &,/*cell iterator for cell*/
      const int &,/*concerning face number in local cell*/
      dealii::FullMatrix<double> &,
      dealii::FullMatrix<double> &,
      dealii::FullMatrix<double> &,
      dealii::FullMatrix<double> &,
      const int &,
      const int &) override {};

 protected:
  using EquationBase<dim>::dat_ptr_;
  using EquationBase<dim>::dofs_per_cell_;
  using EquationBase<dim>::equ_name_;
  using EquationBase<dim>::fv_;
  using EquationBase<dim>::fvf_;
  using EquationBase<dim>::have_reflective_bc_;
  using EquationBase<dim>::is_reflective_bc_;
  using EquationBase<dim>::is_eigen_problem_;
  using EquationBase<dim>::mat_vec_;
  using EquationBase<dim>::n_dir_;
  using EquationBase<dim>::n_group_;
  using EquationBase<dim>::n_q_;
  using EquationBase<dim>::n_qf_;
  using EquationBase<dim>::omega_;
  using EquationBase<dim>::pre_streaming_;
  using EquationBase<dim>::pre_collision_;
  using EquationBase<dim>::scaled_fiss_transfer_;
  using EquationBase<dim>::xsec_;
  void GetGroupCellScalarFlux (std::vector<double> &to_fill, int group);
  //std::map<int, std::unique_ptr<dealii::Vector<double>>> global_angular_flux_;
  std::map<int, dealii::Vector<double>> global_angular_flux_;
};

#endif // BART_SRC_EQUATION_SELF_ADJOINT_ANGULAR_FLUX_H_
