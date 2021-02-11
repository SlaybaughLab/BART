#ifndef BART_SRC_MATERIAL_MATERIAL_PROTOBUF_HPP_
#define BART_SRC_MATERIAL_MATERIAL_PROTOBUF_HPP_

#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/lac/full_matrix.h>

#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>

#include "material.pb.h"
#include "material/material_i.hpp"
#include "utility/utility_functions.hpp"

namespace bart::material {

//! This class reads in and pre-processes material properties.
/*!
 \author Weixiong Zheng
 \date 2017/04~06
 Modified by ablank@berkeley.edu August 2018

 \todo Add functionality to perform eigenvalue decomposition.
 */
class MaterialProtobuf : public MaterialI {
 public:
  using MaterialI::MaterialIDMappedTo;
  /*!
    constructor using map of numerical IDs to Material objects
    do_nda is currently unused, and may be removed in the future
    fissile_ids must be non-empty if is_eigen_problem is given as true
   */
  MaterialProtobuf(const std::unordered_map<int, Material>& materials,
                   bool is_eigen_problem,
                   bool do_nda,
                   int number_of_groups,
                   int number_of_materials);

  MaterialProtobuf(const std::unordered_map<int, std::string>& material_filename_map,
                   bool is_eigen_problem,
                   bool do_nda,
                   int number_of_groups,
                   int number_of_materials);

  /*!
    gets the necessary information from the parameter handler and
    delegates to the other constructor
  */
  explicit MaterialProtobuf(dealii::ParameterHandler& prm);

  //! default destructor
  ~MaterialProtobuf() override = default;

  // Overrides for MaterialI see that file for documentation.
  auto GetFissileIDMap() const -> MaterialIDMappedTo<bool> override;
  auto GetDiffusionCoef() const -> MaterialIDMappedTo<std::vector<double>> override;
  auto GetSigT() const -> MaterialIDMappedTo<std::vector<double>> override;
  auto GetInvSigT() const -> MaterialIDMappedTo<std::vector<double>> override;
  auto GetQ() const -> MaterialIDMappedTo<std::vector<double>> override;
  auto GetQPerSter() const -> MaterialIDMappedTo<std::vector<double>> override;
  auto GetNuSigF() const -> MaterialIDMappedTo<std::vector<double>> override;
  auto GetSigS() const -> MaterialIDMappedTo<dealii::FullMatrix<double>> override;
  auto GetSigSPerSter() const -> MaterialIDMappedTo<dealii::FullMatrix<double>> override;
  auto GetChiNuSigF() const -> MaterialIDMappedTo<dealii::FullMatrix<double>> override;
  auto GetChiNuSigFPerSter() const -> MaterialIDMappedTo<dealii::FullMatrix<double>> override;

  /*!
    asserts that a Material is self-consistent,
    contains the required data,
    and that that data is valid

    The energy groups, sigma_t, and the scattering matrix are always required.
    nu Sigma_f and Chi (normalized to 1) are required for fissile materials.

    If the material contains Q, the Q data will be required to be valid.
    Other properties are not checked except for MultipleDefinition,
    which will not be thrown for UNKNOWN_VECTOR.

    Can throw one of the following,
    with checks for each exception in the order they are listed:
      MultipleDefinition
      MissingProperty
      WrongNumberOfValues
      WrongSign
      EnergyGroupBoundariesNotSorted
      ChiDoesNotSumToOne (only if fission data is required)

    the name parameter is used for identifying the material in the
    generated exception
    by default, CombinedName(material) is used for the name
  */
  static auto CheckValid(const Material& material, bool require_fission_data = false, std::string name = "") -> void;

  /*!
    asserts that the materials in the given map are compatible with one another
    the following things are required to be the same in all of the materials:
      number of groups (throws NumberOfGroupsMismatch if different)
      energy group boundaries  (throws EnergyGroupsMismatch if different)
  */
  static auto CheckConsistent(const MaterialIDMappedTo<Material>& materials) -> void;

 private:
  const MaterialIDMappedTo<Material> materials_; //!< protobuf generated Material object for each material ID

  const bool is_eigen_problem_; //!< Boolean to determine if it's an eigenvalue problem.
  const bool do_nda_; //!< Boolean to determine if nonlinear diffusion acceleration is used.
  const int n_group_; //!< Number of energy groups.
  const int n_material_; //!< Number of materials.
  std::unordered_set<int> fissile_ids_; //!< Set of fissile material IDs.

  MaterialIDMappedTo<bool> is_material_fissile_; //!< Map of material ID to a bool indicating if it is fissile.
  MaterialIDMappedTo<std::vector<double>> diffusion_coef_; //!< Diffusion coefficient of all groups for all materials.
  MaterialIDMappedTo<std::vector<double>> sigt_; //!< \f$\sigma_\mathrm{t}\f$ of all groups for all materials.
  MaterialIDMappedTo<std::vector<double>> inv_sigt_; //!< \f$1/\sigma_\mathrm{t}\f$ of all groups for all materials
  MaterialIDMappedTo<std::vector<double>> chi_; //!< \f$\chi\f$ for all materials.
  MaterialIDMappedTo<std::vector<double>> nusigf_; //!< \f$\nu\sigma_\mathrm{f}\f$ of all groups for all materials.
  MaterialIDMappedTo<std::vector<double>> q_; //!< Fixed source value \f$Q\f$'s of all groups for all materials.
  MaterialIDMappedTo<std::vector<double>> q_per_ster_; //!< Fixed source value \f$Q/(4\pi)\f$'s of all groups for all materials.
  MaterialIDMappedTo<dealii::FullMatrix<double>> sigs_; //!< Scattering transfer matrices for all materials.
  MaterialIDMappedTo<dealii::FullMatrix<double>> sigs_per_ster_; //!< Scattering transfer matrices scaled by \f$4\pi\f$
  /*!
   * \f$\chi\nu\sigma_\mathrm{f}\f$ for all materials. It is pre-computed and transformed into a form of transfer matrix
   * for saving computations. The rows correspond to the incoming energy group and the columns correspond to the outgoing
   * neutron energy group.
   */
  MaterialIDMappedTo<dealii::FullMatrix<double>> chi_nusigf_;
  MaterialIDMappedTo<dealii::FullMatrix<double>> chi_nusigf_per_ster_; //!< \f$\chi\nu\sigma_\mathrm{f}/(4\pi)\f$ for all materials.

  // Validation & population functions
  //! populates is_material_fissile_ with all material IDs and a bool indicating if it is in fissile_ids_
  auto PopulateFissileMap() -> void;
  //! populates all Get-able material data from materials_
  auto PopulateData() -> void;
  //! throws an exception if fissile_ids_ contains an unused material id or if fissile_ids_ is empty in an eigen problem
  auto CheckFissileIDs() const -> void;
  //! throws an exception if n_material_ doesn't equal the number of materials read in
  auto CheckNumberOfMaterials() const -> void;
  //! throws an exception if n_group_ doesn't equal the number of groups in each material
  auto CheckNumberOfGroups() const -> void;
  //! calls checkValid on each material
  auto CheckValidEach() const -> void;
  //! calls CheckConsistent on materials_;
  auto CheckConsistent() const -> void;

  //! Reads file names listed in the parameter handler into an unordered_map<int, std::string>
  static auto ReadMaterialFileNames(dealii::ParameterHandler& prm) -> MaterialIDMappedTo<std::string>;
  /*! parses Materials from the given file names accepts both serialized and human-readable protobuf files
   * consistent with material.proto can throw FailedToFindMaterialFile or FailedToParseMaterialFile */
  static auto ParseMaterials(const MaterialIDMappedTo<std::string>& file_name_map) -> MaterialIDMappedTo<Material>;
  //! returns number, full name, abbreviation, and id of a material as one self-describing string
  static auto CombinedName(const std::pair<int, Material>& id_material_pair,
                           const std::string& delimiter = "|") -> std::string;
  //! returns full name, abbreviation, and id of a material as one self-describing string
  static auto CombinedName(const Material& material, const std::string& delimiter = "|") -> std::string;
  //! extracts all vector properties from the Material can throw MultipleDefinition exception
  static auto GetVectorProperties(const Material& material) ->
  std::unordered_map<Material::VectorId, std::vector<double>, std::hash<int>>;
  /*!  returns the first instance of the specified vector property found in the Material will return an empty vector {}
   * if the property isn't found */
  static auto GetVectorProperty(const Material& material, Material::VectorId property_id) -> std::vector<double>;

  /*! returns the first instance of a SIGMA_S matrix found in the Material will return an empty zero by zero matrix if
   * SIGMA_S isn't found can throw WrongNumberOfValues exception */
  static auto GetScatteringMatrix(const Material& material) -> dealii::FullMatrix<double>;

 public:
  /* dealII macros for declaring exceptions are documented in https://www.dealii.org/9.0.0/doxygen/deal.II/group__Exceptions.html
   * parameters are exception name, argument data types, and output sequence for print_info() */

  DeclExceptionMsg(NoFissileIDs, "At least one material ID must be specified as fissile for eigen problems.");

  DeclException2(FailedToFindMaterialFile, std::string, int, << "Failed to find material file \"" << arg1
      << "\" for material number " << arg2 << " in the current working directory.");

  DeclException2(FailedToParseMaterialFile, std::string, int, << "Failed to parse file \"" << arg1
      << "\" for material number " << arg2 << " as either a human-readable or serialized material file defined by "
                                              "material.proto.");
  
  DeclException2(WrongNumberOfMaterials, unsigned int, int, << "The actual number of materials read in (" << arg1
      << ") does not match the number of materials specified (" << arg2 << ").");

  DeclException3(WrongNumberOfGroups, std::string, unsigned int, int, << "The number_of_groups in material " << arg1
      << " does not match the number of groups in MaterialProtobuf." << " (" << arg2 << " != " << arg3 << ")");

  DeclException2(NumberOfGroupsMismatch, std::string, std::string, << "The number_of_groups in material " << arg1
      << " and material " << arg2 << " are different.");

  DeclException3(EnergyGroupsMismatch, std::string, std::string, unsigned int,
                 << "The energy groups boundaries in material " << arg1 << " and material " << arg2
                 << " differ at index " << arg3 << ".");

  /* invalid material exceptions follow arg1 is always a string identifying the material checked in the order listed
   * here */
  DeclException2(MultipleDefinition, std::string, std::string, << "Material " << arg1 << " has multiple definitions of "
      << arg2 << ", which is not allowed.");

  DeclException2(MissingProperty, std::string, std::string, << "Material " << arg1 << " is missing required property "
      << arg2 << ".");

  DeclException4(WrongNumberOfValues, std::string, std::string, unsigned int, unsigned int,
                 << "In material " << arg1 << ", the number of values under the label " << arg2 << " (" << arg3 << ")"
                     << " is inconsistent with number_of_groups = " << arg4 << ".");

  DeclException4(WrongSign, std::string, std::string, double, unsigned int, << "Material " << arg1
      << " contains a value of the wrong sign in " << arg2 << "." << " (" << arg3 << " at index " << arg4 << ")");

  DeclException4(EnergyGroupBoundariesNotSorted, std::string, double, double, double,
                 << "Energy group boundary vector for material " << arg1 << " contains {" << arg2 << ", " << arg3 << ", "
                     << arg4 << "}, which are out of order.");

  /*! CheckValid expects normalization to machine precision levels, which means experimental Chi data with lower precision
   * would have to be changed before going into the .material file */
  DeclException2(ChiDoesNotSumToOne, std::string, double, << "In material " << arg1
      << ", the sum of values in the Chi vector is " << std::setprecision(std::numeric_limits<double>::digits10 + 1)
      << arg2 << ". It must be normalized to 1.");
};

} // namespace bart::material

#endif // BART_SRC_MATERIAL_MATERIAL_PROTOBUF_HPP_
