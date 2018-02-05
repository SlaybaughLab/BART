#ifndef BART_SRC_COMMON_PROBLEM_DEFINITION_H__
#define BART_SRC_COMMON_PROBLEM_DEFINITION_H__

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/quadrature_lib.h>

#include <iostream>
#include <string>
#include <map>
#include <unordered_map>
#include <vector>
#include <set>


const int k_nmat = 50;
const int k_ngrp = 30;

//! This class performs parameter parsing.
/*!
 This class handles user-defined parameters. For details about parameter parsing,
 refer to <a href="https://www.dealii.org/8.5.0/doxygen/deal.II/classParameterHa
 ndler.html" style="color:blue"><b>ParameterHandler</b></a>.

 \author Weixiong Zheng
 \date 2017/06
 */
class ProblemDefinition
{
public:
  //! Class constructor.
  ProblemDefinition ();

  //! Class destructor.
  ~ProblemDefinition ();

  /*!
   This function process ParameterHandler object using info read from user-provided
   input. Specifically, it declares all possible parameter entries and parse info
   from user-defined input file to prm. After the processing, prm contains all the
   necessary info to define the problem and ready for other classes to retrieve
   the info.

   \param prm ParameterHandler object.
   \return Void.
   */
  static void declare_parameters (dealii::ParameterHandler &prm);
};


#endif  // BART_SRC_COMMON_PROBLEM_DEFINITION_H__
