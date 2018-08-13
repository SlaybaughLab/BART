#ifndef BART_SRC_COMMON_BART_BUILDER_H_
#define BART_SRC_COMMON_BART_BUILDER_H_

#include "../aqdata/aq_base.h"
#include "../mesh/mesh_generator.h"
#include "../material/materials.h"
#include "../equation/equation_base.h"
#include "../iteration/eigen_base.h"
#include "../iteration/mg_base.h"
#include "../iteration/ig_base.h"

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
   equation and NDA if required based on parameters specified in prm.

   \param prm dealii::ParameterHandler object.
   \param fe_ptrs A Hash table containing pointers of FE spaces.
   \return Void.
   */
  template <int dim>
  void BuildFESpaces (const dealii::ParameterHandler &prm,
      std::unordered_map<std::string, dealii::FiniteElement<dim, dim>*>& fe_ptrs);

  //! Function used to build AQ data.
  /*!
   The main functionality is to build angular quadrature data based on parameters
   specified in prm. For 1D, Gauss-Legendre quadrature will be built while in
   multi-D, AQBase will be cast to specific model.

   \param prm dealii::ParameterHandler object.
   \param aq_ptr Angular quadrature pointer that needs to be built.
   \return Void.
   */
  template <int dim>
  void BuildAQ (const dealii::ParameterHandler &prm,
      std::unique_ptr<AQBase<dim>> &aq_ptr);

  //! Function used to build material
  /*!
   The main functionality is to build pointer to object of Materials
   based on parameters specified in prm.

   \param prm dealii::ParameterHandler object.
   \param mat_ptr Materials object pointer.
   \return Void.
   */
  void BuildMaterial (dealii::ParameterHandler &prm,
      std::unique_ptr<Materials> &mat_ptr);

  //! Function used to build material
  /*!
   The same as previous function but returning pointer to material.

   \param prm dealii::ParameterHandler object.
   \return Materials object pointer.
   */
  std::unique_ptr<Materials> BuildMaterial (
      dealii::ParameterHandler &prm);

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

  /*!
   The main functionality is to return a pointer to newly constructed objects of
   EquationBase.

   \param prm dealii::ParameterHandler object.
   \param dat_ptr Resource struct object.
   */
  template <int dim>
  std::unordered_map<std::string, std::unique_ptr<EquationBase<dim>>>BuildEqu (
      const dealii::ParameterHandler &prm,
      std::shared_ptr<FundamentalData<dim>> &dat_ptr);

  template <int dim>
  std::unique_ptr<EigenBase<dim>> BuildEigenItr (
      const dealii::ParameterHandler &prm,
      std::shared_ptr<FundamentalData<dim>> &dat_ptr);

  template <int dim>
  std::unique_ptr<MGBase<dim>> BuildMGItr (const dealii::ParameterHandler &prm,
      std::shared_ptr<FundamentalData<dim>> &dat_ptr);

  template <int dim>
  std::unique_ptr<IGBase<dim>> BuildIGItr (const dealii::ParameterHandler &prm,
      std::shared_ptr<FundamentalData<dim>> &dat_ptr);
};

#endif // BART_SRC_COMMON_BART_BUILDER_H_
