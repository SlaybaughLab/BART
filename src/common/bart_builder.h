#ifndef BART_SRC_COMMON_BART_BUILDER_H_
#define BART_SRC_COMMON_BART_BUILDER_H_

// aq data
#include "../aqdata/aq_base.h"
#include "../aqdata/lsgc.h"

#include <vector>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/fe/fe.h>

//! This class provides builders to build objects in BART.
/*!
 This class provides functionalities to build objects used in BART, e.g. FE
 spaces. The motivation is s.t. builder functions originally living in
 BARTDriver will be separated to increase the clearness. And it will also
 benifit the unit testing.
 \author Weixiong Zheng
 \date 2018/03
 */
template<int dim>
class BARTBuilder {
 public:
  //! Constructor.
  /*!
  SetParams will be invoked to set up parameters.

  \param prm dealii::ParameterHandler object.
  */
  BARTBuilder (dealii::ParameterHandler &prm);

  ~BARTBuilder ();

  //! Function used to initialize all the parameters from user input.
  /*!
  The main functionality is to initialize the parameters after read in user
  inputs.

  \param prm dealii::ParameterHandler object.
  \return Void.
  */
  void SetParams (dealii::ParameterHandler &prm);

  //! Function used to build FE spaces for transport equations.
  /*!
  The main functionality is to produce finite element spaces for transport
  equation and NDA if required.

  \param fe_ptrs A vector containing pointers of FE spaces.
  \return Void.
  */
  void BuildFESpaces (std::vector<dealii::FiniteElement<dim, dim>*> &fe_ptrs);

  //! Function used to build AQ data.
  /*!
  The main functionality is to build angular quadrature data. For 1D, Gauss-Legendre
  quadrature will be built while in multi-D, AQBase will be cast to specific model.

  \param aq_ptr Angular quadrature pointer that needs to be built.
  \return Void.
  */
  void BuildAQ (dealii::ParameterHandler &prm,
      std::unique_ptr<AQBase<dim>> &aq_ptr);

 private:
  bool do_nda_;//!< Boolean to determine if NDA is to be used for accelerations.

  unsigned int p_order_;//!< FE polynomial order.

  int aq_order_;//!< SN order.

  std::string ho_discretization_;//!< Spatial discretization for HO equation.
  std::string nda_discretization_;//!< Spatial discretization for NDA equation.
  std::string aq_name_;//!< Angular quadrature name.
};

#endif // BART_SRC_COMMON_BART_BUILDER_H_
