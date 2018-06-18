#ifndef BART_SRC_EQUATION_EVEN_PARITY_H_
#define BART_SRC_EQUATION_EVEN_PARITY_H_

#include "equation_base.h"

//! This class provides weak formulation for even-parity equation.
/*!
 The governing equation for this class is the even-parity equation. Considering
 one-group steady-state equation with isotropic scattering, we transform the
 angular flux to get a even parity flux:
 \f[
 \psi^+\left(\vec{\Omega}\right)=\frac{\psi\left(\vec{\Omega}\right)+
 \psi\left(-\vec{\Omega}\right)}{2}.
 \f]
 One would realize:
 \f[
 \phi=2\int\limits_{2\pi}d\Omega\ \psi^+.
 \f]
 Skipping all the intermediate steps, one would reach the equation:
 \f[
 \vec{\Omega}\cdot\nabla\frac{1}{\sigma_\mathrm{t}}\vec{\Omega}\cdot\nabla\psi^+
 +
 \sigma_\mathrm{t}\psi^+=
 \frac{\sigma_\mathrm{s}}{4\pi}\phi
 +
 \frac{\nu\sigma_\mathrm{f}}{4\pi k_\mathrm{eff}}\phi\ (\mathrm{or}\ \frac{Q}{4\pi}).
 \f]

 This class implements weak formulation for even-parity equation in both DFEM
 and CFEM. DFEM formulation is based on symmetric interior penalty method
 combined with formulations specifically for degenerate diffusion equation, which
 contains a rank 2 tensor in streaming term. Mathematical details of similar
 formulations are extensively discussed <a href="http://epubs.siam.org/doi/abs/1
 0.1137/S0036142900374111" style="color:blue"><b>here</b></a>. For details about
 even-parity equation, please refer to <a href="http://www.springer.com/cda/cont
 ent/document/cda_downloaddocument/9789048134106-c2.pdf?SGWID=0-0-45-1125054-p17
 3921104"><b>E. E. Lewis's book chapter</b></a>.

 Related documentation could be found in <a href="https://www.dealii.org/8.5.0/
 doxygen/deal.II/" style="color:blue"><b>deal.II documentation</b></a>.

 \author Weixiong Zheng
 \date 2017/05
 */
template<int dim>
class EvenParity : public EquationBase<dim> {
public:
  /*!
   Class constructor.

   \param equation_name An abbreviated name of the equation.
   \param prm ParameterHandler object
   */
  EvenParity (const std::string &equation_name
      const dealii::ParameterHandler &prm,
      std::shared_ptr<FundamentalData<dim>> &dat_ptr);

  /*!
   Class destructor.
   */
  ~EvenParity ();

  /*!
   This function provides pre-assembled matrices for even-parity equation in cell.
   Current implementation does not include the interface parts, that being said,
   only CFEM pre-assembly is supported.

   \return Void.
   */
  void PreassembleCellMatrices ();

  /*!
   This function provides cell-wise integrator for bilinear form. Specifically,
   this overrides EquationBase<dim>'s integrator for even parity equation. Note
   that interface and boundary face terms are not handled in this integrator.

   \param cell Active iterator used to iterating cells.
   \param cell_matrix Local matrix in bilinear form to be modified.
   \param g Group index.
   \param i_dir Direction index.
   \return Void.
   */
  void IntegrateCellBilinearForm (
      typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
      dealii::FullMatrix<double> &cell_matrix,
      const int &g,
      const int &dir);

  /*!
   This function provides cellwise integrator for bilinear form assembly on
   boundary.

   \param cell Active cell iterator containing cell info.
   \param fn Face index in current cell for current boundary face.
   \param cell_matrix Local cell matrix to be modified in current cell.
   \param g Group index.
   \param dir Direction index.
   \return Void.
   */
  void IntegrateBoundaryBilinearForm (
      typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
      const int &fn,/*face number*/
      dealii::FullMatrix<double> &cell_matrix,
      const int &g,
      const int &dir);

  /*!
   This function provides integrator for interface bilinear form assembly in DFEM
   formulations for even parity using interior penalty method.

   \param cell Active cell iterator containing cell info.
   \param neigh Cell iterator for neighboring cell about current face.
   \param fn Face index in current cell for current boundary face.
   \param vi_ui Face matrix from testing interior basis by interior basis.
   \param vi_ue Face matrix from testing exterior basis by interior basis.
   \param ve_ui Face matrix from testing interior basis by exterior basis.
   \param ve_ue Face matrix from testing exterior basis by exterior basis.
   \param g Group index.
   \param dir Direction index.
   \return Void.
   */
  void IntegrateInterfaceBilinearForm (
      typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
      typename dealii::DoFHandler<dim>::cell_iterator &neigh,/*cell iterator for cell*/
      unsigned int &fn,/*concerning face number in local cell*/
      dealii::FullMatrix<double> &vi_ui,
      dealii::FullMatrix<double> &vi_ue,
      dealii::FullMatrix<double> &ve_ui,
      dealii::FullMatrix<double> &ve_ue,
      const int &g,
      const int &dir);

  /*!
   This function provides cellwise integrator for linear form assembly specifically
   for the contribution of scattering.

   \param cell Active cell iterator containing cell info.
   \param cell_rhs Local vector (linear form) to be modified.
   \param sflxes_proc Scalar fluxes for all groups living on current processor.
   \param g Group index.
   \param i_dir Direction index.
   \return Void.
   */
  void IntegrateScatteringLinearForm (
      typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
      dealii::Vector<double> &cell_rhs,
      const int &g,
      const int &dir);

  /*!
   This function provides cellwise integrator for linear form assembly specifically
   for the contribution of fixed source in fixed-source problems or fission in
   eigenvalue problems.

   \param cell Active cell iterator containing cell info.
   \param cell_rhs Local vector (linear form) to be modified.
   \param sflxes_prev Scalar fluxes from previous generation due to fission for
   all groups living on current processor.
   \param g Group index.
   \param i_dir Direction index.
   \return Void.

   \note sflxes_prev will do nothing inside the integrator fixed source problems.
   */
  void IntegrateCellFixedLinearForm (
      typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
      dealii::Vector<double> &cell_rhs,
      const int &g,
      const int &dir);

private:
  //! Polynomial order related part of penalty coefficient.
  double c_penalty_;

  //! Frobenius norms of directional tensors (Rank 2).
  /*!
   Denoting the directions as column vectors with length of dim, the directional
   tensor is defined as:
   \f[
   \mathcal{K}:=\vec{\Omega}\vec{\Omega}^\top,
   \f]
   which is a Rank 2 tensor as well as a symmetric dimxdim matrix.

   For the definition of Fronbenius norm, please refer to <a href="http://mathwo
   rld.wolfram.com/FrobeniusNorm.html" style="color:blue"><b>wolfram page</b></a>.
   */
  std::vector<double> tensor_norms_;
};

#endif // BART_SRC_EQUATION_EVEN_PARITY_H_
