#ifndef BART_SRC_COMMON_PROBLEM_DEFINITION_H_
#define BART_SRC_COMMON_PROBLEM_DEFINITION_H_

#include <deal.II/base/parameter_handler.h>

const int kNMat = 50;
const int kNGrp = 30;

//! This namespace contains functions to perform parameter parsing.
/*!
 This namespace handles user-defined parameters. For details about parameter parsing,
 refer to <a href="https://www.dealii.org/8.5.0/doxygen/deal.II/classParameterHa
 ndler.html" style="color:blue"><b>ParameterHandler</b></a>.

 \author Weixiong Zheng
 \date 2017/06, 2018/03
 */
namespace bparams {
//! ParameterHandler object in the global scope of BART
extern dealii::ParameterHandler GlobPrm;

/*!
 This function process ParameterHandler object using info read from user-provided
 input. Specifically, it declares all possible parameter entries and parse info
 from user-defined input file to ParameterHandler object local_prm in the arguments.
 After the processing, local_prm contains all the necessary info to define the
 problem and ready for other classes to retrieve the info.

 \param local_prm A ParameterHandler object to be used for declaring parameter entries.
 \return Void.
 */
void DeclareParameters (dealii::ParameterHandler &local_prm);

/*!
 This function will call the previous function to declare parameter entries in
 bparams::prm.

 \return Void.
*/
void DeclareParameters ();
}


#endif  // BART_SRC_COMMON_PROBLEM_DEFINITION_H_
