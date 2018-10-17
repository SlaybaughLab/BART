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
   * \f$\varphi\f$, this function integrates the following two terms using the
   * quadrature specified in the problem definition:
   * \f[
   * A(i,j) = 
   * \int_{K}\left(\vec{\Omega}\cdot\nabla\varphi_i\right)\frac{1}{\sigma_t}\left(
   * \vec{\Omega}\cdot\nabla\varphi_j\right) dV +
   * \int_{K}\sigma_t\varphi_i\varphi_j dV
   * \f]
   *
   * where \f$\psi\f$ is the angular flux. Adds the result to position
   * \f$(i,j)\f$ in cell_matrix.
   */
  void IntegrateCellBilinearForm (
      typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
      dealii::FullMatrix<double> &cell_matrix,
      const int &g,
      const int &dir) override;
  
  /*!
   * \brief Integrates the linear scattering terms in the SAAF equation and adds
   *        the values to the vector cell_rhs.
   *           
   * For a given cell in the triangulation, \f$K \in T_K\f$, with basis functions
   * \f$\varphi\f$, this function integrates the following two terms using the
   * quadrature specified in the problem definition:
   * \f[
   * \int_{K}\frac{\sigma_s}{4\pi}\phi\varphi dV + \int_{K}\left(\vec{\Omega}
   * \cdot \nabla \varphi\right)\frac{\sigma_s}{4\pi\sigma_t}\phi dV
   * \f]
   *
   * where \f$\phi\f$ is the scalar flux. Adds the result per cell DOF to the
   * input-output vector cell_rhs.
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
