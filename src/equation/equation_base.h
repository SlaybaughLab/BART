#ifndef __equation_base_h__
#define __equation_base_h__

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/fe/fe_poly.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/conditional_ostream.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <unordered_map>
#include <utility>
#include <vector>

#include "../common/problem_definition.h"
#include "../common/preconditioner_solver.h"
#include "../mesh/mesh_generator.h"
#include "../material/material_properties.h"
#include "../aqdata/aq_base.h"

using namespace dealii;

//! This class provides weak form assembly and physical quantity computation functionalities
/*!
 The governing equation can be presented as
 \f[
 \mathcal{T}\psi=\mathcal{S}\psi+\mathcal{F}\psi\ (\mathrm{or}\
 \frac{Q}{4\pi}),
 \f]
 where \f$\mathcal{T}\f$, \f$\mathcal{S}\f$, \f$\mathcal{F}\f$ and
 \f$Q\f$ are the transport operator (including streaming and collision),
 scattering operator, fission operator and fixed volumetric source, respectively.
 
 
 This class implements abstract functionalities of matrix and vector assembly of
 weak formulation for the governing equation above with finite element method as
 well as physical quantity computations such as fission source and 
 \f$k_\mathrm{eff}\f$. It serves as the base class for any equation involved in 
 BART calculation.
 
 \author Weixiong Zheng
 \date 2017/06~08
 */
template <int dim>
class EquationBase
{
public:
  /*!
   Class constructor.
   
   \param equation_name An abbreviated name of the equation.
   \param msh_ptr An shared_ptr of MeshGenerator<dim> object.
   \param aqd_ptr An shared_ptr of AQBase<dim> object.
   \param mat_ptr An shared_ptr of MaterialProperties object.
   */
  EquationBase (std::string equation_name,
                const ParameterHandler &prm,
                const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
                const std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr,
                const std_cxx11::shared_ptr<MaterialProperties> mat_ptr);
  
  /*!
   Virtual class destructor.
   */
  virtual ~EquationBase ();
  
  /*!
   Virtual function to assemble bilinear form. Inside this function, volumetric,
   boundary and cell interface (if applicable) bilinear forms will be assembled.
   
   \return Void.
   */
  virtual void assemble_bilinear_form ();
  
  /*!
   Virtual function to assemble volumetric and boundary bilinear form. It is not
   recommended to override this function until there's a strong reason. What it
   basically does is per component, we go over all cells on current processor
   and call integrators for volumetric and boundary bilinear form assembly.
   
   \note Component loop needs to be the outer loop to avoid MPI related error
   caused by PETSc.
   
   \return void
   */
  virtual void assemble_volume_boundary_bilinear_form ();
  
  /*!
   Virtual function to assemble interface bilinear form for DFEM. It is not
   recommended to override this function until there's a strong reason. What it
   basically does is per component, we go over all cells on current processor
   and call integrators for cell interface bilinear form assembly.
   
   \note Component loop needs to be the outer loop to avoid MPI related error
   caused by PETSc.
   
   \return void
   */
  virtual void assemble_interface_bilinear_form ();
  
  /*!
   Virtual function to assemble closure terms in NDA weak form.
   
   \param ho_equ_ptr Pointer of HO equation.
   \param do_assembly Boolean to determine if closure terms are assembled. By
   default, it's true.
   \return Void.
   */
  virtual void assemble_closure_bilinear_form
  (std_cxx11::shared_ptr<EquationBase<dim> > ho_equ_ptr,
   bool do_assembly = true);
  
  /*!
   Virtual function to assemble linear forms for a specific group using input 
   flux. Overriding has to be provided. Preassumably, fission source or fixed
   source have been assembled before calling this function.
   
   \param sflxes_proc Scalar fluxes for all groups on current processor.
   \param g Group index.
   \return Void.
   */
  virtual void assemble_linear_form
  (std::vector<Vector<double> > &sflxes_proc,
   unsigned int &g);
  
  /*!
   Virtual function to assemble fixed source or fission source linear forms.
   Overriding has to be provided. This is the step before assembling linear 
   forms.
   
   \return Void.
   */
  virtual void assemble_fixed_linear_form
  (std::vector<Vector<double> > &sflx_prev);
  
  /*!
   Virtual function for preassembling streaming and collision matrices at all
   quadrature points in reference cells. Main intension is to reduce computational
   cost on assembly of system matrices, which can be extremely expensive.
   
   \param cell Active cell iterator containing cell info.
   \param streaming_at_qp Preassembled local streaming matrices of all dirs at
   all quadrature points in reference cell to be modified.
   \param collision_aq_qp Collision (mass) matrices at all quadrature points
   in reference cell to be modified.
   \return Void.
   
   \note If overriding is not provided, the cellwise and boundary integrator has
   to provide ad hoc assembly from scratch.
   */
  virtual void pre_assemble_cell_matrices
  (typename DoFHandler<dim>::active_cell_iterator &cell,
   std::vector<std::vector<FullMatrix<double> > > &streaming_at_qp,
   std::vector<FullMatrix<double> > &collision_at_qp);
  
  /*!
   Virtual function providing integrator for volumetric bilinear form assembly.
   This include both streaming and collision term. Boundary bilinear form, yet,
   is assembled in this integrator.
   
   \param cell Active cell iterator containing cell info.
   \param cell_matrix Local matrix for current cell to be modified.
   \param streaming_at_qp Preassembled local streaming matrices of all dirs at 
   all quadrature points in reference cell.
   \param collision_aq_qp Collision (mass) matrices at all quadrature points
   in reference cell.
   \param g Group index.
   \param i_dir Direction index.
   \return Void.
   
   \note Integrator per call only provides integration for one component of an
   equation specified by group and direction index.
   */
  virtual void integrate_cell_bilinear_form
  (typename DoFHandler<dim>::active_cell_iterator &cell,
   FullMatrix<double> &cell_matrix,
   std::vector<std::vector<FullMatrix<double> > > &streaming_at_qp,
   std::vector<FullMatrix<double> > &collision_at_qp,
   const unsigned int &g,
   const unsigned int &i_dir);
  
  /*!
   Virtual function providing integrator for boundary face bilinear form assembly.
   Overriding has to be provided per equation.
   
   \param cell Active cell iterator containing cell info.
   \param fn Face index in current cell for current face.
   \param cell_matrix Local matrix for current cell to be modified.
   \param g Group index.
   \param i_dir Direction index.
   \return Void.
   */
  virtual void integrate_boundary_bilinear_form
  (typename DoFHandler<dim>::active_cell_iterator &cell,
   unsigned int &fn,/*face number*/
   FullMatrix<double> &cell_matrix,
   const unsigned int &g,
   const unsigned int &i_dir);
  
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
  virtual void integrate_interface_bilinear_form
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
   Virtual function to provide cellwise integrator for linear form assembly
   specifically for the contribution of scattering. It needs to be overriden
   for different derived classes.
   
   \param cell Active cell iterator containing cell info.
   \param cell_rhs Local vector (linear form) to be modified.
   \param sflxes_proc Scalar fluxes for all groups living on current processor.
   \param g Group index.
   \param i_dir Direction index.
   \return Void.
   
   \note Integrator per call only provides integration for one component of an
   equation specified by group and direction index.
   */
  virtual void integrate_scattering_linear_form
  (typename DoFHandler<dim>::active_cell_iterator &cell,
   Vector<double> &cell_rhs,
   std::vector<Vector<double> > &sflx_proc,
   const unsigned int &g,
   const unsigned int &i_dir);
  
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
   \param sflxes_prev Scalar fluxes from previous generation due to fission for
   all groups living on current processor.
   \param g Group index.
   \param i_dir Direction index.
   \return Void.
   
   \note sflxes_prev will do nothing inside the integrator fixed source problems. 
   Integrator per call only provides integration for one component of an 
   equation specified by group and direction index.
   */
  virtual void integrate_cell_fixed_linear_form
  (typename DoFHandler<dim>::active_cell_iterator &cell,
   Vector<double> &cell_rhs,
   std::vector<Vector<double> > &sflx_prev,
   const unsigned int &g,
   const unsigned int &i_dir);
  
  /*!
   Virtual function to provide integrator for linear form assembly on boundary.
   Overriding has to be provided if needed per derived class of EquationBase<dim>.
   
   \param cell Active cell iterator containing cell info.
   \param fn Face index in current cell for current face.
   \param cell_rhses Local vector to be modified for boundary contribution of RHS
   of the equation.
   \param g Group index.
   \param i_dir Direction index.
   \return Void.
   */
  virtual void integrate_boundary_linear_form
  (typename DoFHandler<dim>::active_cell_iterator &cell,
   unsigned int &fn,/*face number*/
   Vector<double> &cell_rhses,
   const unsigned int &g,
   const unsigned int &i_dir);
  
  /*!
   Virtual function to perform one-pass linear solve for all components in target
   group. As diffusion has no concept of direction, this function has to be 
   overriden in such a case. For SN, yet, no overriden is needed unless Krylov
   method is developed.
   
   \param g Group index.
   \return Void.
   */
  virtual void solve_in_group (const unsigned int &g);
  
  /*
  // TODO: if DFEM-NDA is developed, this has to be redesigned
  virtual void prepare_cell_corrections
  (const std::vector<std::vector<Tensor<1, dim> > > &ho_cell_dpsi,
   const std::vector<Tensor<1, dim> > &ho_cell_dphi,
   const std::vector<double> &ho_cell_phi,
   std::vector<Tensor<1, dim> > &cell_corrections);
  
  virtual void prepare_boundary_corrections
  (const std::vector<double> &ho_bd_psi,
   const std::vector<double> &ho_bd_phi,
   std::vector<double> &boundary_corrections);
   */
  
  /*!
   Virtual function to initialize:
   
   (1) Global matrices (PETSc objects) for all components;
   
   (2) Global vectors (PETSc objects), i.e. solutions and right hand side 
   vectors;
   
   (3) \f$\phi\f$ for all groups living on current processor.
   
   \param dsp Sparsity pattern built from BartDriver.
   \param local_dofs Index set for indices of degrees of freedom living on
   current processor.
   \param sflxes_proc \f$\phi\f$ for all groups on current processor.
   */
  virtual void initialize_system_matrices_vectors
  (DynamicSparsityPattern &dsp,
   IndexSet &local_dofs,
   std::vector<Vector<double> > &sflxes_proc);
  
  /*!
   Virtual function to generate moments on current processor and update the old
   moments for all groups. By default, scalar flux will be generated from SN
   angular fluxes.
   
   \param sflxes
   \return Void.
   
   \todo Only scalar fluxes are generated currently. Higher moments generation
   has to be implemented if anisotropic scattering exists.
   */
  virtual void generate_moments
  (std::vector<Vector<double> > &sflxes_proc,
   std::vector<Vector<double> > &sflxes_proc_old);
  
  /*!
   Virtual function to generate scalar flux for a specific group. By default, 
   scalar flux will be generated from SN angular fluxes. \f$\phi\f$ from previous
   iteration will also be updated.
   
   \param sflxes_proc \f$\phi\f$ generated by angular flux.
   \param sflxes_proc_old \f$\phi\f$ from previous iterations to be updated.
   \return Void.
   
   \todo Only scalar flux is generated currently. Moments for the group should be
   implemented for scattering anisotropy.
   */
  virtual void generate_moments
  (Vector<double> &sflx,
   Vector<double> &sflx_old,
   const unsigned int &g);
  
  /*!
   Virtual function to generate scalar fluxes. The intension is s.t. those scalar
   fluxes can be used to generate NDA corrections.
   
   \return Void.
   \note Will only be used when NDA is enabled.
   */
  virtual void generate_moments ();
  
  /*!
   Virtual function to retrieve component index given group and direction indices.
   
   \param direction_index Direction index.
   \param g Group index.
   \return Component index.
   \note Overriding has to be provided for non-SN systems.
   */
  virtual unsigned int get_component_index (unsigned int direction_index,
                                            unsigned int g);
  
  /*!
   Virtual function to retrieve direction index given component index. For
   diffusion-like systems (NDA and diffusion), overriding has to be provided. For
   moment systems, direction could be interpreted as moment index.
   
   \param comp_ind Component index.
   \return Direction index.
   */
  virtual unsigned int get_component_direction (unsigned int comp_ind);
  
  /*!
   Virtual function to retrieve group index given component index. By default, a
   mapping for SN are generated and stored and to be retrieved. For non-SN equations,
   overriding needs to be provided.
   
   \param comp_ind Component index.
   \return Group index.
   */
  virtual unsigned int get_component_group (unsigned int comp_ind);
  
  /*!
   Estimate fission source given scalar fluxes. Given scalar fluxes for all groups
   on current processor, fission source will firstly be calculated on current
   processor and thereafter globally value will be summed up and distributed
   to each processor.
   
   \param sflxes_proc \f$\phi\f$ of all groups on current processor.
   \return Global fission source.
   */
  double estimate_fiss_src (std::vector<Vector<double> > &sflxes_proc);
  
  /*!
   Retrieve cell iterators living on current processor and store them in vector.
   
   \param msh_ptr Pointer of MeshGenerator<dim> object.
   \param dof_handler DoFHandler<dim> object containing all cell info on every
   processor.
   \return Void.
   */
  void initialize_cell_iterators_this_proc
  (const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
   const DoFHandler<dim> &dof_handler);
  
  /*!
   Initialize assembly related objects. These mainly include (but are not limited
   to):
   
   (1) Quadrature rules.
   
   (2) Specific finite element type.
   
   (3) <a><b href="https://www.dealii.org/8.4.0/doxygen/deal.II/classFEValuesBas
   e.html" style="color:blue">FEValuesBase</b></a> related objects.
   
   \param Pointer of polynomial-space finite elements. For details, please refer
   to <a href="https://www.dealii.org/8.5.0/doxygen/deal.II/classFE__Poly.html">
   <b>FE_Poly class page</b></a>.
   \return Void.
   */
  void initialize_assembly_related_objects
  (FE_Poly<TensorProductPolynomials<dim>,dim,dim>* fe);
  
  /*!
   Function to scale \f$\chi\nu\sigma_\mathrm{f}\f$ by \f$k_\mathrm{eff}\f$.
   
   \param keff \f$k_\mathrm{eff}\f$ value.
   \return Void.
   */
  void scale_fiss_transfer_matrices (double keff);
  
  /*!
   Function to retrieve current equation name.
   
   \return A string for the equation name.
   */
  std::string get_equ_name ();
protected:
  /*!
   This function returns the mapping result: (boundary id, direction index)->
   refl. index.
   
   \param boundary_id A integer representing ID of the boundary (ranging from 0 to
   2*dim)
   \param incident_angle_index Incident direction index.
   */
  unsigned int get_reflective_direction_index (unsigned int boundary_id,
                                               unsigned int incident_angle_index);
  
  // "c" in the following quantities means "correction" for NDA use
  //!< Pointer of quadrature rule in cell.
  std_cxx11::shared_ptr<QGauss<dim> > q_rule;
  
  //!< Pointer of quadrature rule on cell face.
  std_cxx11::shared_ptr<QGauss<dim-1> > qf_rule;
  
  //!< Pointer of quadrature rule in cell for NDA correction term.
  std_cxx11::shared_ptr<QGauss<dim> > qc_rule;
  
  //!< Pointer of quadrature rule on cell face for NDA correction term.
  std_cxx11::shared_ptr<QGauss<dim-1> > qfc_rule;
  
  //! Pointer of FEValues object.
  /*!
   In short FEValues and FEFaceValues represent finite element evaluated in
   quadrature points of a cell or on a face. For details, please refer to <a href
   ="https://www.dealii.org/8.5.0/doxygen/deal.II/classFEValuesBase.html" 
   style="color:blue"><b>FEValues page</b></a>.
   */
  std_cxx11::shared_ptr<FEValues<dim> > fv;
  
  //! Pointer of FEFaceValues object.
  std_cxx11::shared_ptr<FEFaceValues<dim> > fvf;
  
  //! Pointer of FEFaceValues object used in DFEM interface bilinear form assembly.
  std_cxx11::shared_ptr<FEFaceValues<dim> > fvf_nei;
  
  //! Pointer of FEValues object for NDA cell correction evaluation and assembly.
  std_cxx11::shared_ptr<FEValues<dim> > fvc;
  
  //! Pointer of FEFaceValues object for NDA face correction evaluation and assembly.
  std_cxx11::shared_ptr<FEFaceValues<dim> > fvfc;
  
  std::string equation_name;//!< String for equation name.
  std::string discretization;//!< String for spatial discretization.
  
  double keff;//!< keff with current generation of neutron.
  double keff_prev_gen;//!< keff from last generation of neutron.
  double fission_source;//!< Fission source with current generation neutron.
  double fission_source_prev_gen;//!< Fission source with previous generation neutron.
  
  bool is_eigen_problem;//!< Boolean to determine if it's eigenvalue problem.
  bool do_nda;//!< Boolean to determine if NDA is performed.
  bool have_reflective_bc;//!< Boolean to determine if problem has reflective BC.
  
  unsigned int n_q;//!< Total number quadrature points in a cell for assembly.
  unsigned int n_qf;//!< Total number of face quadrature points in face assembly.
  unsigned int n_qc;//!< Total number of quadrature points for NDA correction assembly in cells.
  unsigned int n_qfc;//!< Total number of quadrature points for NDA correction assembly on faces.
  unsigned int dofs_per_cell;//!< Total number of degrees of freedom per cell.
  
  unsigned int n_dir;//!< Total number of directions if applicable.
  unsigned int n_azi;//!< Total number of azimuthal angles if applicable.
  unsigned int n_total_vars;//!< Total number of components in current equation.
  unsigned int n_group;//!< Total number of groups.
  unsigned int n_material;//!< Total number of material types.
  unsigned int p_order;//!< Polynomial order for finite elements.
  
  //! Quadrature order for NDA volumetric correction assembly.
  const unsigned int nda_quadrature_order;
  unsigned int global_refinements;//! Number uniform refinements to perform.
  
  std::vector<typename DoFHandler<dim>::active_cell_iterator> local_cells;
  
  //! Local to global indices for all cells.
  std::vector<types::global_dof_index> local_dof_indices;
  
  //! Local to global indices for all cells.
  /*!
   The same as local_dof_indices except this is used for neighboring cells when
   assembling DFEM interface terms.
   */
  std::vector<types::global_dof_index> neigh_dof_indices;
  
  std::vector<Tensor<1, dim> > omega_i;//!< All the directions.
  std::vector<double> wi;//!< All the angular weight.

  //! \f$\sigma_\mathrm{t}\f$ of all groups for all materials.
  std::vector<std::vector<double> > all_sigt;
  
  //! \f$1/\sigma_\mathrm{t}\f$ of all groups for all materials.
  std::vector<std::vector<double> > all_inv_sigt;
  
  //! \f$Q\f$ values of all groups for all materials.
  std::vector<std::vector<double> > all_q;
  
  //! \f$Q/(4\pi)\f$ values of all groups for all materials.
  std::vector<std::vector<double> > all_q_per_ster;
  
  //! \f$\nu\sigma_\mathrm{f}\f$ of all groups for all fissile materials.
  /*!
   \todo Change data type to std::unordered_map.
   */
  std::vector<std::vector<double> > all_nusigf;
  
  //! Scattering matrices for all materials (i.e. \f$\sigma_\mathrm{s,g'\to g}\f$).
  /*!
   \todo Change data type to std::vector<FullMatrix<double> >
   */
  std::vector<std::vector<std::vector<double> > > all_sigs;
  
  //! \f$\sigma_\mathrm{s,g'\to g}/(4\pi)\f$
  std::vector<std::vector<std::vector<double> > > all_sigs_per_ster;
  
  //! \f$\chi\nu\sigma_\mathrm{f}\f$ of all incident and outgoing groups for fissile materials.
  std::vector<std::vector<std::vector<double> > > all_chi_nusigf;
  
  //! \f$\chi\nu\sigma_\mathrm{f}/(4\pi)\f$ for fissile materials.
  std::vector<std::vector<std::vector<double> > > all_chi_nusigf_per_ster;
  
  //! \f$\chi\nu\sigma_\mathrm{f}/k_\mathrm{eff}\f$.
  std::vector<std::vector<std::vector<double> > > scaled_fiss_transfer;
  
  //! \f$\chi\nu\sigma_\mathrm{f}/(4\pi k_\mathrm{eff})\f$.
  std::vector<std::vector<std::vector<double> > > scaled_fiss_transfer_per_ster;
  
  //! Mapping: (group, direction)->component index.
  std::map<std::pair<unsigned int, unsigned int>, unsigned int> component_index;
  
  //! Mapping: (reflective boundary ID, direction)->reflective direction ID.
  std::map<std::pair<unsigned int, unsigned int>, unsigned int> reflective_direction_index;
  
  //! Hash table for mapping: component index->(group index, direction index).
  std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int> > inverse_component_index;
  
  //! Hash table for mapping: boundary ID->if boundary is reflective.
  std::unordered_map<unsigned int, bool> is_reflective_bc;
  
  //! Hash table for mapping: material ID->if material is fissile.
  std::unordered_map<unsigned int, bool> is_material_fissile;
  
  //! ostream on processor with rank==0.
  ConditionalOStream pcout;
private:
  /*!
   Function to process input get necessary parameters for equation assembly.
   Relevant parameters will be retrieved from input pointers and assigned to
   correspoding member variables of this class.
   
   \param msh_ptr Pointer of MeshGenerator<dim> object.
   \param aqd_ptr Pointer of AQBase<dim> object.
   \param mat_ptr Pointer of MaterialProperties object.
   \return Void.
   */
  void process_input (const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
                      const std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr,
                      const std_cxx11::shared_ptr<MaterialProperties> mat_ptr);
  
  //! Pointer of PreconditionerSolver object serving as linear solve facility.
  std_cxx11::shared_ptr<PreconditionerSolver> alg_ptr;
  
  //! System matrices for bilinear forms of all components.
  std::vector<PETScWrappers::MPI::SparseMatrix*> sys_mats;
  
  //! System vectors for linear forms of all components.
  std::vector<PETScWrappers::MPI::Vector*> sys_rhses;
  
  //! System vectors for linear forms for fixed or fission source of all components.
  std::vector<PETScWrappers::MPI::Vector*> sys_fixed_rhses;
  
  //! Solution vectors of all components. They are angular fluxes in SN.
  std::vector<PETScWrappers::MPI::Vector*> sys_aflxes;
  
  //! Solution vectors living on current processor.
  std::vector<Vector<double> > aflxes_proc;
  
  //! HO \f$\phi\f$ of all groups on current processor.
  std::vector<Vector<double> > ho_sflxes_proc;
};

#endif	// define  __equation_base_h__
