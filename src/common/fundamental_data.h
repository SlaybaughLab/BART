#ifndef BART_SRC_COMMON_FUNDAMENTAL_DATA_H_
#define BART_SRC_COMMON_FUNDAMENTAL_DATA_H_

#include "../material/material_properties.h"
#include "../aqdata/aq_base.h"
#include "../mesh/mesh_generator.h"

#include <deal.II/base/parameter_handler.h>

template <int dim>
struct FundamentalData {
  /*!
   Constructor. When constructing, data will be built using builder functions.
   */
  FundamentalData (dealii::ParameterHandler &prm);
  //! Destructor
  ~FundamentalData ();

  std::unique_ptr<AQBase<dim>> aq;//!< Pointer to aq data
  std::unique_ptr<MeshGenerator<dim>> mesh;//!< Pointer to mesh generator
  std::unique_ptr<MaterialProperties> material;//!< Pointer to material properties
};

#endif //BART_SRC_COMMON_FUNDAMENTAL_DATA_H_
