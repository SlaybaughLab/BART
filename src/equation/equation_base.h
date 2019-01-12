#ifndef BART_SRC_EQUATION_EQUATION_BASE_H_
#define BART_SRC_EQUATION_EQUATION_BASE_H_

#include <memory>

#include "../common/computing_data.h"
#include "linear_algebra.h"

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

template <int dim>
class EquationBase {
 public:

  /*!
   * \brief Static factory for classes derived from EquationBase.
   *
   * Instantiates and returns the appropriate equation based
   * on the value specified in the problem as `transport model`.
   */
  static std::unique_ptr<EquationBase<dim>> CreateEquation(
      const dealii::ParameterHandler &prm,
      std::shared_ptr<FundamentalData<dim>> &dat_ptr);

  /*!
   * Enumerator for equations that can be instantiated.
   */  
  enum class EquationType { EvenParity, SAAF };
  
  EquationBase (
      const std::string &equation_name,
      const dealii::ParameterHandler &prm,
      std::shared_ptr<FundamentalData<dim>> &dat_ptr);
  ~EquationBase () = default;

  /*!
   Virtual function to assemble bilinear form. Inside this function, volumetric,
   boundary and cell interface (if applicable) bilinear forms will be assembled.

   \return Void.
   */
  void AssembleBilinearForms();

  /*!
   Virtual function to assemble volumetric and boundary bilinear form. It is not
   recommended to override this function until there's a strong reason. What it
   basically does is per component, we go over all cells on current processor
   and call integrators for volumetric and boundary bilinear form assembly.

   \note Component loop needs to be the outer loop to avoid MPI related error
   caused by PETSc.

   \return void
   */
  virtual void AssembleVolumeBoundaryBilinearForms();

  /*!
   Virtual function providing integrator for volumetric bilinear form assembly.
   This include both streaming and collision term. Boundary bilinear form, yet,
   is assembled in this integrator.

   \param cell Active cell iterator containing cell info.
   \param cell_matrix Local matrix for current cell to be modified.
   \param g Group index.
   \param dir Direction index.
   \return Void.

   \note Integrator per call only provides integration for one component of an
   equation specified by group and direction index.
   */
  virtual void IntegrateCellBilinearForm (
      typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
      dealii::FullMatrix<double> &cell_matrix,
      const int &g,
      const int &dir) = 0;

  /*!
   Virtual function providing integrator for boundary face bilinear form assembly.
   Overriding has to be provided per equation.

   \param cell Active cell iterator containing cell info.
   \param fn Face index in current cell for current face.
   \param cell_matrix Local matrix for current cell to be modified.
   \param g Group index.
   \param dir Direction index.
   \return Void.
   */
  virtual void IntegrateBoundaryBilinearForm (
      typename dealii::DoFHandler<dim>::active_cell_iterator &,
      const int &,
      dealii::FullMatrix<double> &cell_matrix,
      const int &,
      const int &dir) = 0;

  /*!
   Virtual function for preassembling streaming and collision matrices at all
   quadrature points in reference cells. Main intension is to reduce computational
   cost on assembly of system matrices, which can be extremely expensive.

   \return Void.

   \note This function is pure virtual meaning override has to be provided.
   */
  virtual void PreassembleCellMatrices () = 0;

  /*!
   Virtual function to assemble linear forms for a specific group using input
   flux. Overriding has to be provided. Preassumably, fission source or fixed
   source have been assembled before calling this function.

   \param g Group index.
   \return Void.
   */
  virtual void AssembleLinearForms (const int &g);

  /*!
   Virtual function to assemble interface bilinear form for DFEM. What it
   basically does is per component, we go over all cells on current processor
   and call integrators for cell interface bilinear form assembly. Override is
   generally not needed.

   \note Component loop needs to be the outer loop to avoid MPI related error
   caused by PETSc.

   \return void
   */
  virtual void AssembleInterfaceBilinearForms ();

  /*!
   Virtual function providing integrator for interface bilinear form assembly in
   DFEM formulations. Mathematically, this contribute to the numerical flux term
   for DFEM. Generically, we would separate out four terms.

   \param cell Active cell iterator containing cell info.
   \param neigh Cell iterator for neighboring cell about current face.
   \param fn Face index in current cell for current face.
   \param vi_ui Face matrix from testing interior basis by interior basis.
   \param vi_ue Face matrix from testing exterior basis by interior basis.
   \param ve_ui Face matrix from testing interior basis by exterior basis.
   \param ve_ue Face matrix from testing exterior basis by exterior basis.
   \param g Group index.
   \param i_dir Direction index.
   \return Void.

   \note Integrator per call only provides integration for one component of an
   equation specified by group and direction index.
   */
  virtual void IntegrateInterfaceBilinearForm (
      typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
      typename dealii::DoFHandler<dim>::cell_iterator &neigh,/*cell iterator for cell*/
      const int &fn,/*concerning face number in local cell*/
      dealii::FullMatrix<double> &vi_ui,
      dealii::FullMatrix<double> &vi_ue,
      dealii::FullMatrix<double> &ve_ui,
      dealii::FullMatrix<double> &ve_ue,
      const int &g,
      const int &i_dir) = 0;

  /*!
   Virtual function to provide cellwise integrator for linear form assembly
   specifically for the contribution of scattering. It needs to be overriden
   for different derived classes.

   \param cell Active cell iterator containing cell info.
   \param cell_rhs Local vector (linear form) to be modified.
   \param g Group index.
   \param dir Direction index.
   \return Void.

   \note Integrator per call only provides integration for one component of an
   equation specified by group and direction index.
   */
  virtual void IntegrateScatteringLinearForm (
      typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
      dealii::Vector<double> &cell_rhs,
      const int &g,
      const int &dir) = 0;
  /*!
   Virtual function to provide integrator for linear form assembly on boundary.
   Overriding has to be provided if needed per derived class of EquationBase<dim>.

   \param cell Active cell iterator containing cell info.
   \param fn Face index in current cell for current face.
   \param cell_rhs Local vector to be modified for boundary contribution of RHS
   of the equation.
   \param g Group index.
   \param dir Direction index.
   \return Void.
   */
  virtual void IntegrateBoundaryLinearForm (
      typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
      const int &fn,/*face number*/
      dealii::Vector<double> &cell_rhs,
      const int &g,
      const int &dir) = 0;

  /*!
   Virtual function to assemble fixed source or fission source linear forms.
   Overriding has to be provided. This is the step before assembling linear
   forms.

   \return Void.
   */
  virtual void AssembleFixedLinearForms ();

  /*!
   Virtual function to provide cellwise integrator for linear form assembly
   specifically for the contribution of fixed source in fixed-source problems or
   fission in eigenvalue problems. It needs to be overriden for different derived
   classes. It can be represented with generic fission operator \f$\mathcal{F}\f$ as
   \f[
   \left(v,\mathcal{F}(\psi)\right)_\mathcal{D},
   \f]

   where \f$v\f$, the test function, does not necessarily need to belong to the same
   function space as \f$\psi\f$.

   \param cell Active cell iterator containing cell info.
   \param cell_rhs Local vector (linear form) to be modified.
   \param g Group index.
   \param dir Direction index.
   \return Void.

   \note sflxes_prev will do nothing inside the integrator fixed source problems.
   Integrator per call only provides integration for one component of an
   equation specified by group and direction index.
   */
  virtual void IntegrateCellFixedLinearForm (
      typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
      dealii::Vector<double> &cell_rhs,
      const int &g,
      const int &dir) = 0;

  /*!
   Virtual function to perform one-pass linear solve for all components in target
   group. As diffusion has no concept of direction, this function has to be
   overriden in such a case. For SN, yet, no overriden is needed unless Krylov
   method is developed.

   \param g Group index.
   \return Void.
   */
  virtual void SolveInGroup (const int &g);

  virtual void GenerateMoments (
      std::map<std::tuple<int,int,int>, dealii::Vector<double>> &moments,
      std::map<std::tuple<int,int,int>, dealii::Vector<double>> &moments_prev,
      const int &g);

  /*!
   Estimate fission source given scalar fluxes. Given scalar fluxes for all groups
   on current processor, fission source will firstly be calculated on current
   processor and thereafter globally value will be summed up and distributed
   to each processor.

   \return Global fission source.
   */
  double EstimateFissSrc();

  /*!
   Function to scale \f$\chi\nu\sigma_\mathrm{f}\f$ by \f$k_\mathrm{eff}\f$.

   \param keff \f$k_\mathrm{eff}\f$ value.
   \return Void.
   */
  void ScaleFissTransferMatrices(double keff);

  /*!
   Virtual function to retrieve component index given group and direction indices.

   \param g Group index.
   \param dir Direction index.
   \return Component index.
   \note Overriding has to be provided for non-SN systems.
   */
  int GetCompInd (const int &g, const int &dir) const;

  /*!
   Function to retrieve direction index given component index. For
   diffusion-like systems (NDA and diffusion), overriding has to be provided. For
   moment systems, direction could be interpreted as moment index.

   \param comp Component index.
   \return Direction index.
   */
  int GetCompDirInd (const int &comp) const;

  /*!
   Function to retrieve group index given component index. By default, a mapping
   for SN are generated and stored and to be retrieved.

   \param comp Component index.
   \return Group index.
   */
  int GetCompGrpInd (const int &comp) const;

  /*!
   Function to retrieve current equation name.

   \return A string for the equation name.
   */
  std::string GetEquName () const;

  /*!
   Function to retrieve total number of variables of current equation.

   \return Integer the total number of variables.
   */
  int GetNTotalVars () const;

 protected:
  const std::string equ_name_;//!< Name of current equation.

  std::shared_ptr<FundamentalData<dim>> dat_ptr_;//!< Pointer to FundamentalData object
  std::shared_ptr<MatrixVector> mat_vec_;//!< Pointer to MatrixVector object
  std::shared_ptr<XSections> xsec_;//!< Pointer to XSection object

  const bool is_eigen_problem_;//!< Boolean to determine if it's eigenvalue problem.
  const bool do_nda_;//!< Boolean to determine if NDA is performed.
  const bool have_reflective_bc_;//!< Boolean to determine if problem has reflective BC.

  const int p_order_;//!< Polynomial order for current equation.
  const int n_dir_;//!< Total number of directions if applicable.
  const int n_group_;//!< Total number of groups.
  const int n_total_vars_;//!< Total number of components in current equation.

  const std::string discretization_;//!< Discretization of current equation.

  //! Hash table for mapping: boundary ID->if boundary is reflective.
  const std::unordered_map<int, bool> is_reflective_bc_;

  //! Pointer of FEValues object.
  /*!
   In short FEValues and FEFaceValues represent finite element evaluated in
   quadrature points of a cell or on a face. For details, please refer to <a href
   ="https://www.dealii.org/8.5.0/doxygen/deal.II/classFEValuesBase.html"
   style="color:blue"><b>FEValues page</b></a>.
   */
  std::shared_ptr<dealii::FEValues<dim>> fv_;

  //! Pointer of FEFaceValues object.
  std::shared_ptr<dealii::FEFaceValues<dim>> fvf_;

  //! Pointer of FEFaceValues object used in DFEM interface bilinear form assembly.
  std::shared_ptr<dealii::FEFaceValues<dim>> fvf_nei_;

  const int n_q_;
  const int n_qf_;
  const int dofs_per_cell_;//!< Total number of degrees of freedom per cell.

  //! Local to global indices for all cells.
  std::vector<dealii::types::global_dof_index> local_dof_indices_;

  /*!
   The same as local_dof_indices except this is used for neighboring cells when
   assembling DFEM interface terms.
   */
  std::vector<dealii::types::global_dof_index> neigh_dof_indices_;

  //! Preassembled streaming matrices at quadrature points
  std::map<std::pair<int,int>, dealii::FullMatrix<double>>
      pre_streaming_;

  //! Preassembled collision matrices at quadrature points
  std::unordered_map<int, dealii::FullMatrix<double>> pre_collision_;

  //! Hash table containing fission transfer matrix.
  /*!
   \f$\chi\nu\sigma_\mathrm{f}/(4\pi k_\mathrm{eff})\f$ for sn equations or
   \f$\chi\nu\sigma_\mathrm{f}/k_\mathrm{eff}\f$ otherwise
   */
  std::unordered_map<int, dealii::FullMatrix<double>>
      scaled_fiss_transfer_;

  std::vector<dealii::Tensor<1, dim>> omega_;//!< All directions in Tensor<1, dim>
  std::vector<double> w_;//!< All angular weights
  std::map<std::pair<int, int>, int> ref_dir_ind_;

 private:
  /*!
   Function to process input get necessary parameters for equation assembly.
   Relevant parameters will be retrieved from input pointers and assigned to
   correspoding member variables of this class.

   \return Void.
   */
  void ProcessInput ();

  //! Linear algebraic solver/preconditioners object.
  LinearAlgebra lin_alg_;

  //! Mapping: (group, direction)->component index.
  std::map<std::pair<int, int>, int> ho_comp_ind_;

  //! Hash table for mapping: component index->(group index, direction index).
  std::unordered_map<int, std::pair<int,int>> ho_inv_comp_ind_;
};

template <int dim>
std::unordered_map<std::string, std::unique_ptr<EquationBase<dim>>>
GetEquations(const dealii::ParameterHandler &prm,
             std::shared_ptr<FundamentalData<dim>> &dat_ptr);

#endif //BART_SRC_EQUATION_EQUATION_BASE_H_
