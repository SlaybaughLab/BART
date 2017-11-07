#ifndef __problem_definition_h__
#define __problem_definition_h__

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/quadrature_lib.h>

#include <iostream>
#include <string>
#include <map>
#include <unordered_map>
#include <vector>
#include <set>

using namespace dealii;

const unsigned int nmat = 50;
const unsigned int ngrp = 30;
const unsigned int z_levels = 30;
const unsigned int y_levels = 100;

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
  static void declare_parameters (ParameterHandler &prm);
};


#endif  // define  __properties_h__
