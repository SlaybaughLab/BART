#ifndef BART_SRC_COMMON_BART_BUILDER_H_
#define BART_SRC_COMMON_BART_BUILDER_H_

#include <vector>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/fe/fe.h>

template<int dim>
class BartBuilder {
 public:
  BartBuilder (dealii::ParameterHandler &prm);

  ~BartBuilder ();

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

 private:
  bool do_nda_;//!< Boolean to determine if NDA is to be used for accelerations.

  unsigned int p_order_;//!< FE polynomial order.

  std::string ho_discretization_;//!< Spatial discretization for HO equation.
  std::string nda_discretization_;//!< Spatial discretization for NDA equation.
};

#endif // BART_SRC_COMMON_BART_BUILDER_H_
