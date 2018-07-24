#ifndef BART_SRC_COMMON_BART_DRIVER_H_
#define BART_SRC_COMMON_BART_DRIVER_H_

#include "../iteration/iterations.h"
#include "computing_data.h"

template <int dim>
class BARTDriver {
 public:
  BARTDriver (dealii::ParameterHandler &prm);
  ~BARTDriver ();

  /*!
   The function used to run the BART. This is the only interface expected to be
   call from exterior.

   \return Void.
   */
  void DriveBART();

 private:
  /*!
   The function used to make a grid for the problem. It is a wrapper function which
   internally calls MeshGenerator::MakeGrid to make a grid associated with the
   triangulation.

   \return Void.
   */
  void MakeGrid();

  /*!
   Function to initialize:

   (1) Global matrices (PETSc objects) for all components;

   (2) Global vectors (PETSc objects), i.e. solutions and right hand side
   vectors;

   (3) Moments for all groups living on current processor.

   \return Void.
   */
  void InitMatVec();

  /*!
   Function used to perform all the iterative calculations. This function is a
   wrapper function which internally calls Iterations::DoIterations to perform
   necessary calculations.

   \return Void.
   */
  void DoIterations();

  /*!
   This function outputs results. The main functionality is to provide output
   files that can be read by Paraview or Visit.

   Procedure of outputing results includes:
   (1) Output results on current processor to .vtu format files
   (2) Write out a .pvtu files on one processor which has access information for
       all .vtu files.
   In the end, one could either get results for all processors at once by reading
   .pvtu file or access results living on specific processors by reading .vtu files.
   */
  void OutputResults() const;

  int n_proc_;//!< Number of processors used in the calculations
  dealii::Triangulation<dim> tria;//!< Triangulation used in serial calculations

  //! Distributed triangulation for parallel calculations.
  dealii::parallel::distributed::Triangulation<dim> distributed_tria;

  dealii::IndexSet local_owned_dofs_;//!< Index of dofs living on current processor

  //! Index of dofs relevant to current processor
  /*!
   This variable includes not only indices of dofs living on current processor,
   but also those living in ghost cells belonging to current processor. Ghost cells
   could be interpreted as cells living on other processors but neighboring to cells
   in current processor.
   */
  dealii::IndexSet local_relevant_dofs_;

  //! Constraints for hanging nodes.
  /*!
   In radiation transport, constraints are more sophiscated. This one is only for
   setting up system without any real values. remove it in future development
   if adaptive mesh refinements are of interest.
   */
  dealii::ConstraintMatrix dummy_constraints_;

  const int n_group_;//!< Number of groups.
  const bool is_eigen_problem_;//!< Boolean to determine if it's eigenvalue problem.
  const bool do_nda_;//!< Boolean to determine if NDA is performed.
  const std::string ho_equ_name_;//!< Name of the HO equation.
  const std::string ho_discretization_;//!< Spatial discretization for HO.
  const std::string nda_discretization_;//!< Spatial discretization for NDA.
  const std::string output_fname_;//!< Name base for output files.

  //! Pointer to FundamentalData object containing all basic data like mesh and aq.
  std::shared_ptr<FundamentalData<dim>> dat_ptr_;
  Iterations<dim> iter_cls_;//!< Iteration class object
  //! All pointers to equations in Hash table
  std::unordered_map<std::string, std::unique_ptr<EquationBase<dim>>> equ_ptrs_;
};

#endif //BART_SRC_COMMON_BART_DRIVER_H_
