#ifndef BART_SRC_EQUATION_SELF_ADJOINT_ANGULAR_FLUX_H_
#define BART_SRC_EQUATION_SELF_ADJOINT_ANGULAR_FLUX_H_

#include "equation_base.h"

/*!
 * \brief Provides the weak formulation for the Self Adjoint Angular Flux 
 *         (SAAF) equation.
 *
 * This is a second-order formulation of the transport equation, given as:
 * \f[
 * -\vec{\Omega} \cdot \nabla \frac{1}{\Sigma_t}\vec{\Omega}\cdot\nabla\psi
 * + \Sigma_t\psi = \frac{1}{4\pi}\left[\Sigma_s\phi + Q - \vec{\Omega}\cdot
 * \nabla\frac{\Sigma_s\phi + Q}{\Sigma_t}\right]
 * \f]
 * 
 * \author Joshua Rehak
 */

template <int dim>
class SelfAdjointAngularFlux : public EquationBase<dim> {
 public:
  /*! Class constructor */
  SelfAdjointAngularFlux(const std::string equation_name,
                         const dealii::ParameterHandler &prm,
                         std::shared_ptr<FundamentalData<dim>> &data_ptr);
  ~SelfAdjointAngularFlux() = default;


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
   * where \f$\psi\f$ is the angular flux. Adds the result to position
   * \f$(i,j)\f$ in cell_matrix
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

  // NOT DONE YET
  
  

  void IntegrateBoundaryBilinearForm (
      typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
      const int &fn,
      dealii::FullMatrix<double> &cell_matrix,
      const int &g,
      const int &dir) override {};

  void IntegrateInterfaceBilinearForm (
      typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
      typename dealii::DoFHandler<dim>::cell_iterator &neigh,/*cell iterator for cell*/
      const int &fn,/*concerning face number in local cell*/
      dealii::FullMatrix<double> &vi_ui,
      dealii::FullMatrix<double> &vi_ue,
      dealii::FullMatrix<double> &ve_ui,
      dealii::FullMatrix<double> &ve_ue,
      const int &g,
      const int &i_dir) override {};
  void IntegrateBoundaryLinearForm (
      typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
      const int &fn,/*face number*/
      dealii::Vector<double> &cell_rhs,
      const int &g,
      const int &dir) override {};  

 protected:
  using EquationBase<dim>::xsec_;
  using EquationBase<dim>::pre_streaming_;
  using EquationBase<dim>::pre_collision_;
  using EquationBase<dim>::omega_;
  using EquationBase<dim>::equ_name_;
  using EquationBase<dim>::fv_;
  using EquationBase<dim>::dat_ptr_;
  using EquationBase<dim>::mat_vec_;
  using EquationBase<dim>::n_group_;
  using EquationBase<dim>::n_q_;
  using EquationBase<dim>::n_dir_;
  using EquationBase<dim>::dofs_per_cell_;
  
  dealii::FullMatrix<double> CellCollisionMatrix (int q);
  dealii::FullMatrix<double> CellStreamingMatrix (int q, int dir);
};

#endif // BART_SRC_EQUATION_SELF_ADJOINT_ANGULAR_FLUX_H_
