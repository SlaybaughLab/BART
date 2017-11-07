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
 This class implements abstract functionalities of matrix and vector assembly as
 well as physical quantity computations such as fission source and keff. It serves
 as the base class for any equation involved in BART calculation.
 
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
   Class constructor. The difference from previous constructor is that no AQBase
   pointer is needed. It was designed for diffusion-like system.
   
   \param equation_name An abbreviated name of the equation.
   \param msh_ptr An shared_ptr of MeshGenerator<dim> object.
   \param mat_ptr An shared_ptr of MaterialProperties object.
   
   \todo Discuss the necessity of this constructor. We actually can use the first
   constructor simply by passing empty AQBase pointer for diffusion.
   */
  EquationBase (std::string equation_name,
                const ParameterHandler &prm,
                const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
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
  
  virtual void assemble_volume_boundary_bilinear_form ();
  
  virtual void assemble_interface_bilinear_form ();
  
  void assemble_closure_bilinear_form
  (std_cxx11::shared_ptr<EquationBase<dim> > ho_equ_ptr,
   bool do_assembly = true);
  
  virtual void assemble_linear_form
  (std::vector<Vector<double> > &sflx_this_proc,
   unsigned int &g);
  
  virtual void assemble_fixed_linear_form
  (std::vector<Vector<double> > &sflx_prev);
  
  virtual void pre_assemble_cell_matrices
  (typename DoFHandler<dim>::active_cell_iterator &cell,
   std::vector<std::vector<FullMatrix<double> > > &streaming_at_qp,
   std::vector<FullMatrix<double> > &collision_at_qp);
  
  virtual void integrate_cell_bilinear_form
  (typename DoFHandler<dim>::active_cell_iterator &cell,
   FullMatrix<double> &cell_matrix,
   std::vector<std::vector<FullMatrix<double> > > &streaming_at_qp,
   std::vector<FullMatrix<double> > &collision_at_qp,
   const unsigned int &g,
   const unsigned int &i_dir);
  
  virtual void integrate_boundary_bilinear_form
  (typename DoFHandler<dim>::active_cell_iterator &cell,
   unsigned int &fn,/*face number*/
   FullMatrix<double> &cell_matrix,
   const unsigned int &g,
   const unsigned int &i_dir);
  
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
  
  virtual void integrate_scattering_linear_form
  (typename DoFHandler<dim>::active_cell_iterator &cell,
   Vector<double> &cell_rhs,
   std::vector<Vector<double> > &sflx_proc,
   const unsigned int &g,
   const unsigned int &i_dir);
  
  virtual void integrate_cell_fixed_linear_form
  (typename DoFHandler<dim>::active_cell_iterator &cell,
   Vector<double> &cell_rhs,
   std::vector<Vector<double> > &sflx_prev,
   const unsigned int &g,
   const unsigned int &i_dir);
  
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
  
  virtual void initialize_system_matrices_vectors
  (DynamicSparsityPattern &dsp,
   IndexSet &local_dofs,
   std::vector<Vector<double> > &sflxes_proc);
  
  virtual void generate_moments
  (std::vector<Vector<double> > &sflxes,
   std::vector<Vector<double> > &sflxes_old);
  
  virtual void generate_moments
  (Vector<double> &sflx,
   Vector<double> &sflx_old,
   const unsigned int &g);
  
  virtual void generate_moments ();
  
  // override this three functions in derived classes
  // these functions have to be redefined when using diffusion, NDA, PN, SPN
  virtual unsigned int get_component_index (unsigned int incident_angle_index, unsigned int g);
  virtual unsigned int get_component_direction (unsigned int comp_ind);
  virtual unsigned int get_component_group (unsigned int comp_ind);
  
  double estimate_fiss_src (std::vector<Vector<double> > &sflxes_proc);
  
  void initialize_cell_iterators_this_proc
  (const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
   const DoFHandler<dim> &dof_handler);
  
  void initialize_assembly_related_objects
  (FE_Poly<TensorProductPolynomials<dim>,dim,dim>* fe);
  
  void scale_fiss_transfer_matrices (double keff);
  
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
  std_cxx11::shared_ptr<QGauss<dim> > q_rule;//!< Pointer to quadrature rule in cell.
  std_cxx11::shared_ptr<QGauss<dim-1> > qf_rule;
  std_cxx11::shared_ptr<QGauss<dim> > qc_rule;
  std_cxx11::shared_ptr<QGauss<dim-1> > qfc_rule;
  std_cxx11::shared_ptr<FEValues<dim> > fv;
  std_cxx11::shared_ptr<FEFaceValues<dim> > fvf;
  std_cxx11::shared_ptr<FEFaceValues<dim> > fvf_nei;
  std_cxx11::shared_ptr<FEValues<dim> > fvc;
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
  std::vector<types::global_dof_index> local_dof_indices;
  std::vector<types::global_dof_index> neigh_dof_indices;
  
  std::vector<Tensor<1, dim> > omega_i;//!< All the directions.
  std::vector<double> wi;//!< All the angular weight.
  std::vector<std::vector<double> > all_sigt;
  std::vector<std::vector<double> > all_inv_sigt;
  std::vector<std::vector<double> > all_q;
  std::vector<std::vector<double> > all_q_per_ster;
  std::vector<std::vector<double> > all_nusigf;
  std::vector<std::vector<std::vector<double> > > all_sigs;
  std::vector<std::vector<std::vector<double> > > all_sigs_per_ster;
  std::vector<std::vector<std::vector<double> > > all_chi_nusigf;
  std::vector<std::vector<std::vector<double> > > all_chi_nusigf_per_ster;
  std::vector<std::vector<std::vector<double> > > scaled_fiss_transfer_per_ster;
  std::vector<std::vector<std::vector<double> > > scat_scaled_fiss_transfer_per_ster;
  std::vector<std::vector<std::vector<double> > > scaled_fiss_transfer;
  
  std::map<std::pair<unsigned int, unsigned int>, unsigned int> component_index;
  std::map<std::pair<unsigned int, unsigned int>, unsigned int> reflective_direction_index;
  std::map<std::vector<unsigned int>, unsigned int> relative_position_to_id;
  std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int> > inverse_component_index;
  std::unordered_map<unsigned int, bool> is_reflective_bc;
  std::unordered_map<unsigned int, bool> is_material_fissile;
  
  std::set<unsigned int> fissile_ids;
  
  ConditionalOStream pcout;
private:
  void setup_system ();
  void generate_globally_refined_grid ();
  void report_system ();
  void print_angular_quad ();
  
  // void setup_lo_system();
  void setup_boundary_ids ();
  void get_cell_mfps (unsigned int &material_id, double &cell_dimension,
                      std::vector<double> &local_mfps);
  
  void process_input (const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
                      const std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr,
                      const std_cxx11::shared_ptr<MaterialProperties> mat_ptr);
  
  std_cxx11::shared_ptr<PreconditionerSolver> alg_ptr;
  
  // related objects for current equation: matrices, vectors
  std::vector<PETScWrappers::MPI::SparseMatrix*> sys_mats;
  std::vector<PETScWrappers::MPI::Vector*> sys_rhses;
  std::vector<PETScWrappers::MPI::Vector*> sys_fixed_rhses;
  std::vector<PETScWrappers::MPI::Vector*> sys_aflxes;
  std::vector<Vector<double> > aflxes_proc;
  std::vector<Vector<double> > ho_sflxes_proc;
};

#endif	// define  __equation_base_h__
