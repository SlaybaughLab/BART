#ifndef BART_SRC_ITERATION_ITERATION_BASE_H_
#define BART_SRC_ITERATION_ITERATION_BASE_H_

#include <deal.II/base/parameter_handler.h>

#include "../equation/equation_base.h"

//! This class provides common functionalities for iteration related classes.
/*!
 Currently, it is mainly for providing functionalities on estimating vector
 differences.

 \author Weixiong
 \date 2017/08~10
 \todo Implement iteration counting for every type of iteration in calculations.
 */
template <int dim>
class IterationBase {
 public:
  /*!
   Class constructor.

   \param prm Const dealii::ParameterHandler object.
   */
  IterationBase (const dealii::ParameterHandler &prm,
      std::shared_ptr<FundamentalData<dim>> &dat_ptr);

  //! Virtual class destructor.
  virtual ~IterationBase () = default;

 protected:
  /**
   * Function to measure the relative difference between two sets of PETSc
   * Vectors
   *
   * \param Two sets of PETSc Vectors with the same length
   * \return relative difference in vector l1 norm measure
   */
  double EstimatePhiDiff (
      std::vector<dealii::PETScWrappers::MPI::Vector*> &phis1,
      std::vector<dealii::PETScWrappers::MPI::Vector*> &phis2);

  /** Function to measure the relative difference between two PETSc
   * Vectors
   *
   * \param Two PETSc Vectors with the same length
   * \return relative difference in vector l1 norm measure
   */
  double EstimatePhiDiff (
      dealii::PETScWrappers::MPI::Vector* phi1,
      dealii::PETScWrappers::MPI::Vector* phi2);

  /**
   * Function to measure the relative difference between two PETSc MPI Vectors.
   *
   * \param Two sets of PETSc Vectors with the same length
   * \return relative difference in vector l1 norm measure
   */
  double EstimatePhiDiff (
      std::map<std::tuple<int,int,int>, dealii::Vector<double>> &phis1,
      std::map<std::tuple<int,int,int>, dealii::Vector<double>> &phis2);

  /**
   * Function to measure the relative difference between two sets of deal.II
   * Vectors. Note that though all vectors live on current processor, broadcast
   * is performed to obtain global results.
   *
   * \param Two deal.II Vectors with the same length
   * \return relative difference in vector l1 norm measure
   */
  double EstimatePhiDiff (
      dealii::Vector<double> &phi1, dealii::Vector<double> &phi2);

  const int n_group_;//!< Number of groups.
  const bool is_eigen_problem_;//!< Boolean to determine if it's eigenvalue problem.
  const bool do_nda_;//!< Boolean to determine if NDA is used.

  double total_calculation_time_; /**< total time for calculations+assemblies*/
  int ct_ho_iters_; /**< HO iteration counts*/
  int ct_nda_iters_; /**< NDA iteration counts*/

  std::string ho_equ_name_;//!< High-order equation name
  std::string iter_equ_name_;//!< Name of equation used in iterations

  //! moments from previou iteration
  std::map<std::tuple<int,int,int>, dealii::Vector<double>> moments_prev_;
  std::shared_ptr<FundamentalData<dim>> dat_ptr_;
  std::shared_ptr<MatrixVector> mat_vec_;//!< Pointer to MatrixVector object
};

#endif // BART_SRC_ITERATION_ITERATION_BASE_H_
