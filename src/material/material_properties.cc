#include <deal.II/base/exceptions.h>

#include "material_properties.h"

MaterialProperties::MaterialProperties(std::unordered_map<int, std::string> materials_as_strings,
                                        bool is_eigen_problem,
                                        bool do_nda,
                                        int number_of_groups,
                                        int number_of_materials,
                                        std::unordered_set<int> fissile_ids /* = {} */)
    : is_eigen_problem_(is_eigen_problem),
      do_nda_(do_nda),
      n_group_(number_of_groups),
      n_material_(number_of_materials),
      fissile_ids_(fissile_ids) {
    CheckFissileIDs();
    PopulateFissileMap(); //required for CheckValidEach()

    ParseMaterials(materials_as_strings);

    CheckNumberOfMaterials();
    CheckValidEach();
    CheckNumberOfGroups();
    CheckConsistent();

    PopulateData();
}

MaterialProperties::MaterialProperties(dealii::ParameterHandler& prm)
    : MaterialProperties(
        ReadMaterialStrings(prm),
        prm.get_bool("do eigenvalue calculations"),
        prm.get_bool("do NDA"),
        prm.get_integer("number of groups"),
        prm.get_integer("number of materials"),
        ReadFissileIDs(prm)) {

}

MaterialProperties::~MaterialProperties() {}

///!TODO implement this
std::unordered_map<int, std::string> MaterialProperties::ReadMaterialStrings(dealii::ParameterHandler& prm)
{
  //prm.leave_subsection()
  return {};
}

///!TODO implement this
std::unordered_set<int> MaterialProperties::ReadFissileIDs(dealii::ParameterHandler& prm)
{
  //prm.leave_subsection()
  return {};
}

void MaterialProperties::ParseMaterials(const std::unordered_map<int, std::string>& materials_as_strings)
{
  ///!TODO
  // #include <google/protobuf/text_format.h>; TextFormat::ParseFromString
  // try both human-readable and serialized, throw FailedToParseMaterialString
}

void MaterialProperties::PopulateFissileMap()
{
  ///!TODO
  // do not subtract 1 like in dev branch implementation
}

void MaterialProperties::PopulateData()
{
  ///!TODO implement this
  /* consistent with the dev branch implementation,
   * unordered maps for fixed source data will be left empty for an eigen problem and
   * materials not labeled fissile will be excluded from the unordered_maps for fission data
  */
  //Q assumed to be zero if not given in the Material
}

void MaterialProperties::CheckFissileIDs() const
{
  AssertThrow(!(is_eigen_problem_ && fissile_ids_.empty()), NoFissileIDs());

  for (const int& id : fissile_ids_)
  {
    AssertThrow(materials_.count(id) > 0, FissileIDInvalid(id));
  }
}

void MaterialProperties::CheckNumberOfMaterials() const
{
  AssertThrow(n_material_ > 0 && (unsigned int)(n_material_) == materials_.size(),
    WrongNumberOfMaterials(n_material_, materials_.size()));
}

void MaterialProperties::CheckNumberOfGroups() const
{
  for (const std::pair<int, Material>& mat_pair : materials_)
  {
    const Material& mat = mat_pair.second;
    AssertThrow(n_group_ > 0 && (unsigned int)(n_group_) == mat.number_of_groups(),
      WrongNumberOfGroups(mat_pair.first, CombinedName(mat), mat.number_of_groups(), n_group_));
  }
}

void MaterialProperties::CheckValid(const Material& material, bool require_fission_data /* = false*/)
{
  ///!TODO implement this, make sure error message contains material id and specifically what's wrong
  /*!
    the energy groups, sigma_t, and the scattering matrix are always required
    nu Sigma_f and Chi (normalized to 1) are required for fissile materials
  */
  /*
    check for:
      MissingProperty
      MultipleDefinition
      WrongNumberOfValues
      ChiDoesNotSumToOne
      EnergyGroupBoundariesNotSorted
      WrongSign
  */
}

void MaterialProperties::CheckValidEach() const
{
  for (const std::pair<int, Material>& mat_pair : materials_)
  {
    try
    {
      CheckValid(mat_pair.second, is_material_fissile_.at(mat_pair.first));
    }
    catch (const dealii::ExceptionBase& e)
    {
      throw ExceptionWithID(e, mat_pair.first);
    }
  }
}

void MaterialProperties::CheckConsistent(const std::unordered_map<int, Material>& materials)
{
  ///!TODO implement this
  /*
    will maybe print a warning without failing if multiple IDs map to one material
      is printing a warning kosher?
    will throw an exception if two materials have differrent energy group boundaries
  */

  // throw NumberOfGroupsMismatch, EnergyGroupsMismatch
}

void MaterialProperties::CheckConsistent() const
{
  CheckConsistent(materials_);
}

std::string MaterialProperties::CombinedName(const Material& material, const std::string& delimiter /* = "|" */)
{
  const std::string q = "\"";
  return q + material.full_name() + q + delimiter + q + material.abbreviation() + q + delimiter + q + material.id() + q;
}
