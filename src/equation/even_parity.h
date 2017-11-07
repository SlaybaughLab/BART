#ifndef __even_parity__
#define __even_parity__

#include "equation_base.h"

//! This class provides weak formulation for even-parity equation.
/*!
 This class implements weak formulation for even-parity equation in both DFEM
 and CFEM. DFEM formulation is based on symmetric interior penalty method 
 combined with formulations specifically for degenerate diffusion equation, which
 contains a rank 2 tensor in streaming term. Mathematical details of similar
 formulations are extensively discussed <a href="http://epubs.siam.org/doi/abs/1
 0.1137/S0036142900374111" style="color:blue"><b>here</b></a>.
 
 Related documentation could be found in <a href="https://www.dealii.org/8.5.0/
 doxygen/deal.II/" style="color:blue"><b>deal.II documentation</b></a>
 
 \author Weixiong Zheng
 \date 2017/05
 */
template<int dim>
class EvenParity : public EquationBase<dim>
{
public:
  /*!
   Class constructor.
   
   \param equation_name An abbreviated name of the equation.
   \param msh_ptr An shared_ptr of MeshGenerator<dim> object.
   \param aqd_ptr An shared_ptr of AQBase<dim> object.
   \param mat_ptr An shared_ptr of MaterialProperties object.
   */
  EvenParity (std::string equation_name,
              const ParameterHandler &prm,
              const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
              const std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr,
              const std_cxx11::shared_ptr<MaterialProperties> mat_ptr);
  
  /*!
   Class destructor.
   */
  ~EvenParity ();
  
  /*!
   This function provides pre-assembled matrices for even-parity equation in cell.
   Current implementation does not include the interface parts, that being said,
   only CFEM pre-assembly is supported.
   
   \param cell Active iterator used to iterating cells.
   \param streaming_at_qp Vector of streaming matrices per quadrature point.
   \param collision_at_qp Vector of mass matrices per quadrature point.
   \return Void.
   */
  void pre_assemble_cell_matrices
  (typename DoFHandler<dim>::active_cell_iterator &cell,
   std::vector<std::vector<FullMatrix<double> > > &streaming_at_qp,
   std::vector<FullMatrix<double> > &collision_at_qp);
  
  /*!
   This function provides cell-wise integrator for bilinear form. Specifically,
   this overrides EquationBase<dim>'s integrator for even parity equation. Note
   that interface and boundary face terms are not handled in this integrator.
   
   \param cell Active iterator used to iterating cells.
   \param cell_matrix Local matrix in bilinear form to be modified.
   \param streaming_at_qp Preassembled streaming matrix per quadrature point.
   \param collision_at_qp Preassembled collision matrix per quadrature point.
   \param g Group index.
   \param i_dir Direction index.
   \return Void.
   */
  void integrate_cell_bilinear_form
  (typename DoFHandler<dim>::active_cell_iterator &cell,
   FullMatrix<double> &cell_matrix,
   std::vector<std::vector<FullMatrix<double> > > &streaming_at_qp,
   std::vector<FullMatrix<double> > &collision_at_qp,
   const unsigned int &g,
   const unsigned int &i_dir);
  
  /*!
   This function provides cellwise integrator for bilinear form assembly on
   boundary.
   
   \param cell Active cell iterator containing cell info.
   \param fn Face index in current cell for current boundary face.
   \param cell_matrix Local cell matrix to be modified in current cell.
   \param g Group index.
   \param i_dir Direction index.
   \return Void.
   */
  void integrate_boundary_bilinear_form
  (typename DoFHandler<dim>::active_cell_iterator &cell,
   unsigned int &fn,/*face number*/
   FullMatrix<double> &cell_matrix,
   const unsigned int &g,
   const unsigned int &i_dir);
  
  /*!
   This function provides integrator for interface bilinear form assembly in DFEM
   formulations.
   
   \param cell Active cell iterator containing cell info.
   \param neigh Cell iterator for neighboring cell about current face.
   \param fn Face index in current cell for current boundary face.
   \param vi_ui Face matrix from testing interior basis by interior basis.
   \param vi_ue Face matrix from testing exterior basis by interior basis.
   \param ve_ui Face matrix from testing interior basis by exterior basis.
   \param ve_ue Face matrix from testing exterior basis by exterior basis.
   \param g Group index.
   \param i_dir Direction index.
   \return Void.
   */
  void integrate_interface_bilinear_form
  (typename DoFHandler<dim>::active_cell_iterator &cell,
   typename DoFHandler<dim>::cell_iterator &neigh,/*cell iterator for cell*/
   unsigned int &fn,/*concerning face number in local cell*/
   FullMatrix<double> &vi_ui,
   FullMatrix<double> &vi_ue,
   FullMatrix<double> &ve_ui,
   FullMatrix<double> &ve_ue,
   const unsigned int &g,
   const unsigned int &i_dir);
  
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
  void integrate_scattering_linear_form
  (typename DoFHandler<dim>::active_cell_iterator &cell,
   Vector<double> &cell_rhs,
   std::vector<Vector<double> > &sflxes_proc,
   const unsigned int &g,
   const unsigned int &i_dir);

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
  void integrate_cell_fixed_linear_form
  (typename DoFHandler<dim>::active_cell_iterator &cell,
   Vector<double> &cell_rhs,
   std::vector<Vector<double> > &sflxes_prev,
   const unsigned int &g,
   const unsigned int &i_dir);
  
private:
  //! Polynomial order related part of penalty coefficient.
  double c_penalty;
  
  //! Frobenius norms of directional tensors (Rank 2).
  std::vector<double> tensor_norms;
};

#endif // __even_parity__
