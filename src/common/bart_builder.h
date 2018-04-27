#ifndef BART_SRC_COMMON_BART_BUILDER_H_
#define BART_SRC_COMMON_BART_BUILDER_H_

// aq data
#include "../aqdata/aq_base.h"
#include "../aqdata/lsgc.h"

#include <vector>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/fe/fe.h>

//! This namespace provides builders to build objects in BART.
/*!
 This namespace provides functionalities to build objects used in BART, e.g. FE
 spaces. The motivation is s.t. builder functions originally living in
 BARTDriver will be separated to increase the clearness. And it will also
 benifit the unit testing.
 \author Weixiong Zheng
 \date 2018/03
 */
namespace bbuilders {
  //! Function used to build FE spaces for transport equations.
  /*!
  The main functionality is to produce finite element spaces for transport
  equation and NDA if required.

  \param fe_ptrs A vector containing pointers of FE spaces.
  \return Void.
  */
  template <int dim>
  void BuildFESpaces (const dealii::ParameterHandler &prm,
      std::vector<dealii::FiniteElement<dim, dim>*> &fe_ptrs);

  //! Function used to build AQ data.
  /*!
  The main functionality is to build angular quadrature data. For 1D, Gauss-Legendre
  quadrature will be built while in multi-D, AQBase will be cast to specific model.

  \param aq_ptr Angular quadrature pointer that needs to be built.
  \return Void.
  */
  template <int dim>
  void BuildAQ (const dealii::ParameterHandler &prm,
      std::unique_ptr<AQBase<dim>> &aq_ptr);
};

#endif // BART_SRC_COMMON_BART_BUILDER_H_
