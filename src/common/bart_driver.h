#ifndef __bart_driver_h__
#define __bart_driver_h__

#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/fe/fe_poly.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/numerics/data_out.h>

#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../mesh/mesh_generator.h"
#include "../material/material_properties.h"
#include "../aqdata/aq_base.h"
#include "../equation/equation_base.h"
#include "../iteration/mg_base.h"
#include "../iteration/ig_base.h"
#include "../iteration/eigen_base.h"
#include "../iteration/iterations.h"

using namespace dealii;

/*!
 This class operate BART. The functionalities include:
 
 (1) Build necessary components for BART, such as AQ data, mesh, equations etc. For
     This functionality, please see documentation for BartDriver<dim>::build_basis.
 
 (2) Setting up system.
 
 (3) Reporting system on screen.
 
 (4) Operation calculations at highest level: instantiating Iterations and perform
     iterations.
 
 (5) Output results for visualization.
 
 \author Weixiong Zheng
 \date 2017/08
 */
template <int dim>
class BartDriver
{
public:
  /*!
   Class constructor.
   
   \param prm ParameterHandler object containing all user defined parameters.
   */
  BartDriver (ParameterHandler &prm);
  ~BartDriver ();//!< Destructor
  
  //! Public function used to run all the functions inside it.
  /*!
   This is the only function callable outside of the class.
   */
  void run ();
private:
  /*!
   A function invoked in BartDriver constructor to build basis of BART. The main
   functionalities includes:
   (1) Building pointer to Iteration class.
   (2) Downcasting AQBase and make angular quadrature.
   (3) Building pointer to MeshGenerator class for generating mesh.
   (4) Building pointer to MaterialProperties class for generating material properties
   (5) Downcasting finite element type
   (6) Downcasting EquationBase
   (7) Downcasting InGroupBase, MGBase, EigenBase
   
   \param prm ParameterHandler object containing all parameters.
   \return Void.
   */
  void build_basis (ParameterHandler &prm);
  
  /*!
   This function initializes system settings. Main functionalities include:
   (1) Initializing dof_handler, local_dofs and relevant_dofs.
   (2) Initializing, producing and distributing sparsity pattern
   (3) Initializing cell iterators on current iterators, assembly related objects
       and system matrices and vectors (PETSc related objects)
   */
  void setup_system ();
  
  //! Function to print features such as SN order, number of groups etc. on screen.
  void report_system ();
  
  /*!
   This function outputs results. The main functionality is to provide output 
   files that can be read by Paraview or Visit.
   
   Procedure of outputing results includes:
   (1) Output results on current processor to .vtu format files
   (2) Write out a .pvtu files on one processor which has access information for
       all .vtu files.
   */
  void output_results () const;
  
  /*!
   Function builds pointer of EquationBase object, i.e. instance of space-angle 
   solver. The main functionality is to perform downcasting from base class to 
   derived class based on derived class type defined in ParameterHandler object. 
   For robustness reason, shared_ptr instead of raw pointer is used in the 
   casting process.
   
   \param equation_name A string defining the name of the desired equation.
   \param prm A ParameterHandler object containing all parameters processed from
   user input.
   \return shared_ptr pointing to EquationBase instance.
   */
  std_cxx11::shared_ptr<EquationBase<dim> > build_equation
  (std::string equation_name,
   const ParameterHandler &prm,
   const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
   const std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr,
   const std_cxx11::shared_ptr<MaterialProperties> mat_ptr);
  
  //! Function used to build pointer to instance of AQBase's derived class
  /*!
   The main functionality is to perform downcasting from base class to derived
   class based on derived class type defined in ParameterHandler object. For
   robustness reason, shared_ptr instead of raw pointer is used in the casting
   process.
   
   \param prm A ParameterHandler object containing all parameters processed from
   user input.
   \return shared_ptr pointing to AQBase instance.
   */
  std_cxx11::shared_ptr<AQBase<dim> > build_aq_model (ParameterHandler &prm);
  
  //! Function used to build pointer to instance of EigenBase's derived class
  /*!
   The main functionality is to perform downcasting from base class to derived
   class based on derived class type defined in ParameterHandler object. For 
   robustness reason, shared_ptr instead of raw pointer is used in the casting 
   process.
   
   \param prm A ParameterHandler object containing all parameters processed from
   user input.
   \return shared_ptr pointing to EigenBase instance.
   */
  std_cxx11::shared_ptr<EigenBase<dim> > build_eigen_iterations (const ParameterHandler &prm);
  
  //! Function used to build pointer to instance of MGBase's derived class
  std_cxx11::shared_ptr<MGBase<dim> > build_mg_iterations (const ParameterHandler &prm);
  
  //! Function used to build pointer to instance of InGroupBase's derived class
  std_cxx11::shared_ptr<IGBase<dim> > build_ig_iterations (const ParameterHandler &prm);
  
  //! Finite element type.
  /*!
   The finite element type will be specified by downcasting to either FE_DGQ or
   FE_Q.
   
   \todo It would be necessary to modify it to vectors of fe pointers if CMFD is
   of interest.
   */
  FE_Poly<TensorProductPolynomials<dim>,dim,dim>* fe;
  
  parallel::distributed::Triangulation<dim> triangulation;//!< Triangulation in distrubted system
  
  //! DoFHandler instance
  /*!
   Could be interpreted as a container containing smart pointers to accessors of
   cells (cell iterators) and relevant cell properties such as material_id, if cell
   is at boundary and dof mapping between local cell and global system etc.
   */
  DoFHandler<dim> dof_handler;
  
  IndexSet local_dofs;//!< Index of dofs living on current processor
  
  //! Index of dofs relevant to current processor
  /*!
   This variable includes not only indices of dofs living on current processor,
   but also those living in ghost cells belonging to current processor. Ghost cells
   could be interpreted as cells living on other processors but neighboring to cells
   in current processor.
   */
  IndexSet relevant_dofs;
  
  //! Constraints for hanging nodes.
  /*!
   This will be useful if adaptive local refinement is of interest. For now,
   this variable is of no practical use.
   */
  ConstraintMatrix constraints;
  ConditionalOStream pcout;//!< Conditional ostream, i.e. ostream on one processor.
  
  std_cxx11::shared_ptr<Iterations<dim> > itr_ptr;//!< Pointer of Iterations instance.
  std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr;//!< Pointer of MeshGenerator instance.
  std_cxx11::shared_ptr<MaterialProperties> mat_ptr;//!< Pointer of MaterialProperties instance.
  std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr;//!< Pointer of AQBase instance.
  std_cxx11::shared_ptr<EigenBase<dim> > eig_ptr;//!< Pointer of EigenBase instance.
  std_cxx11::shared_ptr<MGBase<dim> > mg_ptr;//!< Pointer of MGBase instance.
  std_cxx11::shared_ptr<IGBase<dim> > ig_ptr;//!< Pointer of InGroupBase instance.
  
  //! A vector of pointers of EquationBase instances.
  /*!
   The rule is if NDA is not used, this vector only has one component for transport
   model. Otherwise, std::vector::front() contains the pointer of HO equation and
   std::vector::back() contains the pointer of NDA equation.
   */
  std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > equ_ptrs;
  
  std::string transport_model;//!< Name of transport model.
  std::string ho_linear_solver_name;//!< Algebraic solver name.
  std::string ho_preconditioner_name;//!< Algebraic preconditioner name.
  std::string discretization;//!< Spatial discretization method.
  std::string namebase;//< Output filename base.
  
  double keff;//!< keff result in eigenvalue problems
  
  bool is_eigen_problem;//!< Boolean to determine if a problem is about eigenvalue.
  bool do_nda;//!< Boolean to determine if NDA is to be used for accelerations.
  bool have_reflective_bc;//!< Boolean to determine if there is any BC are reflective.
  
  unsigned int n_dir;//!< Total number of directions in SN calculations.
  unsigned int n_azi;//!< SN order.
  unsigned int n_total_ho_vars;//!< Total number of components.
  unsigned int n_group;//!< Total number of energy groups.
  unsigned int n_material;//!< Total number of material types.
  unsigned int p_order;//!< Polynomial order.
  unsigned int global_refinements;//!< Total number of global refinements.
  
  /*!
   The scalar fluxes for all groups on current processor. Values will be assigned
   in the process of calculations. It will also be the output scalar fluxes in the 
   member method output_results
   */
  std::vector<Vector<double> > sflxes_proc;
};

#endif	// define  __bart_driver_h__
