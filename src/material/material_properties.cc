#include "material_properties.h"

MaterialProperties::MaterialProperties (dealii::ParameterHandler &prm)
    : is_eigen_problem_(prm.get_bool("do eigenvalue calculations")),
      do_nda_(prm.get_bool("do NDA")),
      n_group_(prm.get_integer("number of groups")),
      n_material_(prm.get_integer("number of materials")) {
  ProcessMaterialProperties(prm);
}

MaterialProperties::~MaterialProperties () {}

void MaterialProperties::ProcessMaterialProperties(
    dealii::ParameterHandler &prm) {      
}
