#ifndef BART_SRC_MATERIAL_MATERIAL_PROPERTIES_H_
#define BART_SRC_MATERIAL_MATERIAL_PROPERTIES_H_

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/lac/full_matrix.h>

#include <google/protobuf/text_format.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

#include "../../material.pb.h"

#include "../common/numbers.h"

//! This class read in and pre-process material properties.
/*!
 \author Weixiong Zheng
 \date 2017/04~06

 \todo Change scattering data from std::vector<std::vector<std::vector<double> > >
 to std::vector<Vector<double> > for scattering matrix eigenvalue decomposition.
 \todo Add functionality to perform eigenvalue decomposition.
 */
class MaterialProperties {
 public:
  /*!
    constructor using map of numerical IDs to Material objects
   */
  MaterialProperties(const std::unordered_map<int, Material>& materials,
                      bool is_eigen_problem,
                      bool do_nda /*!< currently unused, may be removed in the future */,
                      int number_of_groups,
                      int number_of_materials,
                      const std::unordered_set<int>& fissile_ids = {} /*!< only used in eigen problems */);

  /*!
    gets the necessary information from the parameter handler and delegates to the other constructor
  */
  MaterialProperties(dealii::ParameterHandler& prm);

  //! Class destructor.
  ~MaterialProperties();

  //! returns an unordered_map from material ID to a boolean for whether it is fissile
  std::unordered_map<int, bool> GetFissileIDMap() const;

  //! A function to retrieve all \f$\sigma_\mathrm{t}\f$ for all groups.
  std::unordered_map<int, std::vector<double>> GetSigT() const;

  //! A function to retrieve all \f$1/\sigma_\mathrm{t}\f$ for all groups.
  std::unordered_map<int, std::vector<double>> GetInvSigT() const;

  //! A function to retrieve all fixed source value \f$Q\f$'s for all groups.
  std::unordered_map<int, std::vector<double>> GetQ() const;

  //! A function to retrieve all \f$Q/(4\pi)\f$'s for all groups.
  std::unordered_map<int, std::vector<double>> GetQPerSter() const;

  //! A function to retrieve all \f$\nu\sigma_\mathrm{f}\f$'s.
  std::unordered_map<int, std::vector<double>> GetNuSigF() const;

  //! A function to retrieve all scattering transfer matrices.
  std::unordered_map<int, dealii::FullMatrix<double>> GetSigS() const;

  //! A function to retrieve all scattering transfer matrices scaled by \f$4\pi\f$.
  std::unordered_map<int, dealii::FullMatrix<double>> GetSigSPerSter() const;

  //! A function to retrieve \f$\chi\nu\sigma_\mathrm{f}\f$ for all fissile materials.
  std::unordered_map<int, dealii::FullMatrix<double>> GetChiNuSigF() const;

  //! A function to retrieve \f$\chi\nu\sigma_\mathrm{f}/(4\pi)\f$ for all fissile materials.
  std::unordered_map<int, dealii::FullMatrix<double>> GetChiNuSigFPerSter() const;

  /*!
    tests if a Material is self-consistent, contains the required data,
    and if that data is valid

    The energy groups, sigma_t, and the scattering matrix are always required.
    nu Sigma_f and Chi (normalized to 1) are required for fissile materials.

    If the material contains Q, the Q data will be required to be valid.
    Other properties are not checked except for MultipleDefinition,
    which will not be thrown for UNKNOWN_VECTOR.

    Can throw one of the following, with checks for each exception in the order they are listed:
      MultipleDefinition
      MissingProperty
      WrongNumberOfValues
      WrongSign
      EnergyGroupBoundariesNotSorted
      ChiDoesNotSumToOne (only if fission data is required)

    name parameter is used for identifying the material in the generated exception
    by default, CombinedName(material) is used for the name
  */
  static void CheckValid(const Material& material, bool require_fission_data = false, std::string name = "");

  /*!
    tests if the given map holds materials that are all compatible with one another
    the following things are required to be the same in all of the materials:
      number of groups
      energy group boundaries
  */
  static void CheckConsistent(const std::unordered_map<int, Material>& materials);

 private:
  //! protobuf generated Material object for each material ID
  const std::unordered_map<int, Material> materials_;
  const bool is_eigen_problem_;//!< Boolean to determine if it's eigenvalue problem.
  const bool do_nda_;//!< Boolean to determine if NDA is used.
  const int n_group_;//!< Number of groups.
  const int n_material_;//!< Number of materials.

  //! Hash table to determine if a material is fissile.
  std::unordered_map<int, bool> is_material_fissile_;

  std::unordered_set<int> fissile_ids_;//!< Set of fissile material IDs.

  //! \f$\sigma_\mathrm{t}\f$ of all groups for all materials.
  std::unordered_map<int, std::vector<double>> sigt_;

  //! \f$1/\sigma_\mathrm{t}\f$ of all groups for all materials
  std::unordered_map<int, std::vector<double>> inv_sigt_;

  /*!
   \f$\chi\f$ for all materials. A mistake when designing this is that it was
   treated as other properties. Yet, it is group independent.

   \todo Change this to std::vector<double> chi, or using Hash table.
   */
  std::unordered_map<int, std::vector<double>> chi_;

  //! \f$\nu\sigma_\mathrm{f}\f$ of all groups for all materials.
  std::unordered_map<int, std::vector<double>> nusigf_;

  //! Fixed source value \f$Q\f$'s of all groups for all materials.
  std::unordered_map<int, std::vector<double>> q_;

  //! Fixed source value \f$Q/(4\pi)\f$'s of all groups for all materials.
  std::unordered_map<int, std::vector<double>> q_per_ster_;

  //! Scattering transfer matrices for all materials.
  std::unordered_map<int, dealii::FullMatrix<double>> sigs_;

  //! Scattering transfer matrices scaled by \f$4\pi\f$
  std::unordered_map<int, dealii::FullMatrix<double>> sigs_per_ster_;

  //! \f$\chi\nu\sigma_\mathrm{f}\f$ for all materials.
  /*!
   It is pre-computed and transformed into a form of transfer matrix for saving
   computations.
   */
  std::unordered_map<int, dealii::FullMatrix<double>> chi_nusigf_;

  //! \f$\chi\nu\sigma_\mathrm{f}/(4\pi)\f$ for all materials.
  std::unordered_map<int, dealii::FullMatrix<double>> chi_nusigf_per_ster_;

  //! populates is_material_fissile_ with all material IDs and a bool indicating if it is in fissile_ids_
  void PopulateFissileMap();

  //! populates all Get-able material data from materials_
  void PopulateData();

  //! throws an exception if fissile_ids_ contains an unused material id or if fissile_ids_ is empty in an eigen problem
  void CheckFissileIDs() const;

  //! throws an exception if n_material_ doesn't equal the number of materials read in
  void CheckNumberOfMaterials() const;

  //! throws an exception if n_group_ doesn't equal the number of groups in each material
  void CheckNumberOfGroups() const;

  //! calls checkValid on each material
  void CheckValidEach() const;

  //! calls CheckConsistent on materials_;
  void CheckConsistent() const;

  //! Reads file names listed in the parameter handler into an unordered_map<int, std::string>
  static std::unordered_map<int, std::string> ReadMaterialFileNames(dealii::ParameterHandler& prm);

  //! Reads fissile material ids listed in the parameter handler
  static std::unordered_set<int> ReadFissileIDs(dealii::ParameterHandler& prm);

  /*!
    parses Materials from the given file names
    accepts both serialized and human-readable protobuf files consistent with material.proto
    can throw FailedToFindMaterialFile or FailedToParseMaterialFile
  */
  static std::unordered_map<int, Material> ParseMaterials(const std::unordered_map<int, std::string>& file_name_map);

  //! returns number, full name, abbreviation, and id of a material as one self-describing string
  static std::string CombinedName(const std::pair<int, Material>& id_material_pair, const std::string& delimiter = "|");

  //! returns full name, abbreviation, and id of a material as one self-describing string
  static std::string CombinedName(const Material& material, const std::string& delimiter = "|");

  //! extracts all vector properties from the Material, can throw MultipleDefinition exception
  static std::unordered_map<Material::VectorId, std::vector<double>> GetVectorProperties(const Material& material);

  /*! 
    returns the first instance of the specified vector property found in the Material
    will return an empty vector {} if the property isn't found
  */
  static std::vector<double> GetVectorProperty(const Material& material, Material::VectorId property_id);

  /*! 
    returns the first instance of a SIGMA_S matrix found in the Material
    will return an empty zero by zero matrix if SIGMA_S isn't found
    can throw WrongNumberOfValues exception
  */
  static dealii::FullMatrix<double> GetScatteringMatrix(const Material& material);

  //! returns the sum of the given values calculated using an algorithm that minimizes rounding error
  static double PreciseSum(const std::vector<double>& values);

 public:
  /*
    dealII macros for declaring exceptions are documented in https://www.dealii.org/9.0.0/doxygen/deal.II/group__Exceptions.html
    parameters are exception name, argument data types, and output sequence for print_info()
  */

  DeclExceptionMsg(NoFissileIDs, "At least one material ID must be specified as fissile for eigen problems.");
  
  DeclException1(FissileIDInvalid, int,
    << "Material ID " << arg1 << " was specified as fissile, but no material with ID " << arg1 << " exists.");

  DeclException2(FailedToFindMaterialFile, std::string, int,
    << "Failed to find material file \"" << arg1 << "\" for material number " << arg2
    << " in the current working directory.");

  DeclException2(FailedToParseMaterialFile, std::string, int,
    << "Failed to parse file \"" << arg1 << "\" for material number " << arg2
    << " as either a human-readable or serialized material file defined by material.proto.");
  
  DeclException2(WrongNumberOfMaterials, unsigned int, int,
    << "The actual number of materials read in (" << arg1
    << ") does not match the number of materials specified (" << arg2 << ").");

  DeclException3(WrongNumberOfGroups, std::string, unsigned int, int,
    << "The number_of_groups in material " << arg1
    << " does not match the number of groups in MaterialProperties."
    << " (" << arg2 << " != " << arg3 << ")");

  DeclException2(NumberOfGroupsMismatch, std::string, std::string,
    << "The number_of_groups in material " << arg1 << " and material " << arg2 << " are different.");

  DeclException3(EnergyGroupsMismatch, std::string, std::string, unsigned int,
    << "The energy groups boundaries in material " << arg1 << " and material " << arg2 << " differ at index " << arg3 << ".");

  // invalid material exceptions follow; arg1 is always a string identifying the material; checked in the order listed here

  DeclException2(MultipleDefinition, std::string, std::string,
    << "Material " << arg1 << " has multiple definitions of " << arg2 << ", which is not allowed.");

  DeclException2(MissingProperty, std::string, std::string,
    << "Material " << arg1 << " is missing required property " << arg2 << ".");

  DeclException4(WrongNumberOfValues, std::string, std::string, unsigned int, unsigned int,
    << "In material " << arg1 << ", the number of values under the label " << arg2 << " (" << arg3 << ")"
    << " is inconsistent with number_of_groups = " << arg4 << ".");

  DeclException4(WrongSign, std::string, std::string, double, unsigned int,
    << "Material " << arg1 << " contains a value of the wrong sign in " << arg2 << "."
    << " (" << arg3 << " at index " << arg4 << ")");

  DeclException4(EnergyGroupBoundariesNotSorted, std::string, double, double, double,
    << "Energy group boundary vector for material " << arg1
    << " contains {" << arg2 << ", " << arg3 << ", " << arg4 << "}, which are out of order.");

  /*!
    CheckValid expects normalization to machine precision levels, which means experimental Chi data with
    lower precision would have to be changed before going into the .material file
  */
  DeclException2(ChiDoesNotSumToOne, std::string, double,
    << "In material " << arg1 << ", the sum of values in the Chi vector is "
    << std::setprecision(std::numeric_limits<double>::digits10 + 1) << arg2 << ". It must be normalized to 1.");
};

#endif //BART_SRC_MATERIAL_MATERIAL_PROPERTIES_H_
