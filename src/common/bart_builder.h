#ifndef BART_SRC_COMMON_BART_BUILDER_H_
#define BART_SRC_COMMON_BART_BUILDER_H_

#include "../mesh/mesh_generator.h"

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

  //! Function used to build mesh
  /*!
   The main functionality is to build pointer to object of MeshGenerator<dim>
   based on parameters specified in prm.

   \param prm dealii::ParameterHandler object.
   \param msh_ptr MeshGenerator object pointer.
   \return Void.
   */
  template <int dim>
  void BuildMesh (dealii::ParameterHandler &prm,
      std::unique_ptr<MeshGenerator<dim>> &msh_ptr);

  /*!
   The same as previous function but returning mesh pointer.

   \param prm dealii::ParameterHandler object.
   \return MeshGenerator object pointer.
   */
  template <int dim>
  std::unique_ptr<MeshGenerator<dim>> BuildMesh (
      dealii::ParameterHandler &prm);

}

#endif // BART_SRC_COMMON_BART_BUILDER_H_
