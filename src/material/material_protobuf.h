#ifndef BART_SRC_MATERIAL_MATERIAL_PROTOBUF_H_
#define BART_SRC_MATERIAL_MATERIAL_PROTOBUF_H_

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
#include "../utility/utility_functions.h"

#include "../common/numbers.h"
#include "material_base.h"

//! This class reads in and pre-processes material properties.
/*!
 \author Weixiong Zheng
 \date 2017/04~06
 Modified by ablank@berkeley.edu August 2018

 \todo Add functionality to perform eigenvalue decomposition.
 */

class MaterialProtobuf : public MaterialBase {
 public:
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

  /*!
    returns an unordered_map from material ID to a
    boolean that is true if the material was labeled fissile
  */
  std::unordered_map<int, bool> GetFissileIDMap() const override;

  std::unordered_map<int, std::vector<double>> GetDiffusionCoef() const override;
  
  //! Returns all \f$\sigma_\mathrm{t}\f$ for all groups.
  std::unordered_map<int, std::vector<double>> GetSigT() const override;

  //! Returns all \f$1/\sigma_\mathrm{t}\f$ for all groups.
  std::unordered_map<int, std::vector<double>> GetInvSigT() const override;

  //! Returns all fixed source value \f$Q\f$'s for all groups.
  std::unordered_map<int, std::vector<double>> GetQ() const override;

  //! Returns all \f$Q/(4\pi)\f$'s for all groups.
  std::unordered_map<int, std::vector<double>> GetQPerSter() const override;

  //! Returns all \f$\nu\sigma_\mathrm{f}\f$'s.
  std::unordered_map<int, std::vector<double>> GetNuSigF() const override;

  //! Returns all scattering transfer matrices.
  std::unordered_map<int, dealii::FullMatrix<double>>
  GetSigS() const override;

  //! Returns all scattering transfer matrices scaled by \f$4\pi\f$.
  std::unordered_map<int, dealii::FullMatrix<double>>
  GetSigSPerSter() const override;

  //! Returns fission transfer matrix \f$\chi\nu\sigma_\mathrm{f}\f$ for all fissile materials.
  /*! Returns the fission transfer matrix for all fissile materials. Entries
      are given by:
      \f[
      \mathbf{F} = \begin{bmatrix}
      \nu_0\sigma_{f,0}\chi_0 & \nu_0\sigma_{f,0}\chi_1 & \cdots & \nu_0\sigma_{f,0}\chi_G \\
      \nu_1\sigma_{f,1}\chi_0 & \nu_1\sigma_{f,1}\chi_1 & \cdots & \nu_1\sigma_{f,1}\chi_G \\
      \vdots & \vdots & \ddots & \vdots \\
      \nu_G\sigma_{f,G}\chi_0 & \nu_G\sigma_{f,G}\chi_1 & \cdots & \nu_G\sigma_{f,G}\chi_G
      \end{bmatrix}
      \f]
  */
  std::unordered_map<int, dealii::FullMatrix<double>>
  GetChiNuSigF() const override;

  //! Returns \f$\chi\nu\sigma_\mathrm{f}/(4\pi)\f$ for all fissile materials.
  std::unordered_map<int, dealii::FullMatrix<double>>
  GetChiNuSigFPerSter() const override;

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
  static void CheckValid(const Material& material,
                         bool require_fission_data = false,
                         std::string name = "");

  /*!
    asserts that the materials in the given map are compatible with one another
    the following things are required to be the same in all of the materials:
      number of groups (throws NumberOfGroupsMismatch if different)
      energy group boundaries  (throws EnergyGroupsMismatch if different)
  */
  static void CheckConsistent(const std::unordered_map<int, Material>& materials);

 private:
  //! protobuf generated Material object for each material ID
  const std::unordered_map<int, Material> materials_;

  //! Boolean to determine if it's an eigenvalue problem.
  const bool is_eigen_problem_;

  //! Boolean to determine if nonlinear diffusion acceleration is used.
  const bool do_nda_;

  //! Number of energy groups.
  const int n_group_;

  //! Number of materials.
  const int n_material_;

  //! Set of fissile material IDs.
  std::unordered_set<int> fissile_ids_;

  //! Map of material ID to a bool indicating if it is fissile.
  std::unordered_map<int, bool> is_material_fissile_;

  //! Diffusion coefficient of all groups for all materials.
  std::unordered_map<int, std::vector<double>> diffusion_coef_;
  
  //! \f$\sigma_\mathrm{t}\f$ of all groups for all materials.
  std::unordered_map<int, std::vector<double>> sigt_;

  //! \f$1/\sigma_\mathrm{t}\f$ of all groups for all materials
  std::unordered_map<int, std::vector<double>> inv_sigt_;

  //! \f$\chi\f$ for all materials.
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

  /*!
    \f$\chi\nu\sigma_\mathrm{f}\f$ for all materials.
    It is pre-computed and transformed into a form of transfer matrix
    for saving computations.
    The rows correspond to the incoming energy group and
    the columns correspond to the outgoing neutron energy group.
   */
  std::unordered_map<int, dealii::FullMatrix<double>> chi_nusigf_;

  //! \f$\chi\nu\sigma_\mathrm{f}/(4\pi)\f$ for all materials.
  std::unordered_map<int, dealii::FullMatrix<double>> chi_nusigf_per_ster_;

  /*!
    populates is_material_fissile_ with all material IDs and
    a bool indicating if it is in fissile_ids_
  */
  void PopulateFissileMap();

  //! populates all Get-able material data from materials_
  void PopulateData();

  /*!
    throws an exception if fissile_ids_ contains an unused material id or
    if fissile_ids_ is empty in an eigen problem
  */
  void CheckFissileIDs() const;

  /*!
    throws an exception if n_material_ doesn't equal
    the number of materials read in
  */
  void CheckNumberOfMaterials() const;

  /*!
    throws an exception if n_group_ doesn't equal
    the number of groups in each material
  */
  void CheckNumberOfGroups() const;

  //! calls checkValid on each material
  void CheckValidEach() const;

  //! calls CheckConsistent on materials_;
  void CheckConsistent() const;

  /*!
    Reads file names listed in the parameter handler into an
    unordered_map<int, std::string>
  */
  static std::unordered_map<int, std::string>
  ReadMaterialFileNames(dealii::ParameterHandler& prm);

  /*!
    parses Materials from the given file names
    accepts both serialized and human-readable protobuf files
    consistent with material.proto
    can throw FailedToFindMaterialFile or FailedToParseMaterialFile
  */
  static std::unordered_map<int, Material>
  ParseMaterials(const std::unordered_map<int, std::string>& file_name_map);

  /*!
    returns number, full name, abbreviation, and id of a material
    as one self-describing string
  */
  static std::string
  CombinedName(const std::pair<int, Material>& id_material_pair,
               const std::string& delimiter = "|");

  /*!
    returns full name, abbreviation, and id of a material
    as one self-describing string
  */
  static std::string
  CombinedName(const Material& material,
               const std::string& delimiter = "|");

  /*!
    extracts all vector properties from the Material
    can throw MultipleDefinition exception
  */
  static std::unordered_map<Material::VectorId, std::vector<double>, std::hash<int>>
  GetVectorProperties(const Material& material);

  /*! 
    returns the first instance of the specified vector property
    found in the Material
    will return an empty vector {} if the property isn't found
  */
  static std::vector<double>
  GetVectorProperty(const Material& material, Material::VectorId property_id);

  /*! 
    returns the first instance of a SIGMA_S matrix found in the Material
    will return an empty zero by zero matrix if SIGMA_S isn't found
    can throw WrongNumberOfValues exception
  */
  static dealii::FullMatrix<double>
  GetScatteringMatrix(const Material& material);

 public:
  /*
    dealII macros for declaring exceptions are documented in
    https://www.dealii.org/9.0.0/doxygen/deal.II/group__Exceptions.html
    parameters are exception name,
    argument data types,
    and output sequence for print_info()
  */

  DeclExceptionMsg(NoFissileIDs,
    "At least one material ID must be specified as fissile for eigen problems.");

  DeclException2(FailedToFindMaterialFile,
    std::string, int,
    << "Failed to find material file \"" << arg1 << "\" for material number "
    << arg2 << " in the current working directory.");

  DeclException2(FailedToParseMaterialFile,
    std::string, int,
    << "Failed to parse file \"" << arg1 << "\" for material number " << arg2
    << " as either a human-readable or serialized material file defined by material.proto.");
  
  DeclException2(WrongNumberOfMaterials,
    unsigned int, int,
    << "The actual number of materials read in (" << arg1
    << ") does not match the number of materials specified (" << arg2 << ").");

  DeclException3(WrongNumberOfGroups,
    std::string, unsigned int, int,
    << "The number_of_groups in material " << arg1
    << " does not match the number of groups in MaterialProtobuf."
    << " (" << arg2 << " != " << arg3 << ")");

  DeclException2(NumberOfGroupsMismatch,
    std::string, std::string,
    << "The number_of_groups in material " << arg1
    << " and material " << arg2 << " are different.");

  DeclException3(EnergyGroupsMismatch,
    std::string, std::string, unsigned int,
    << "The energy groups boundaries in material " << arg1 <<
    " and material " << arg2 << " differ at index " << arg3 << ".");

  /*
    invalid material exceptions follow
    arg1 is always a string identifying the material
    checked in the order listed here
  */

  DeclException2(MultipleDefinition,
    std::string, std::string,
    << "Material " << arg1 << " has multiple definitions of "
    << arg2 << ", which is not allowed.");

  DeclException2(MissingProperty,
    std::string, std::string,
    << "Material " << arg1 << " is missing required property "
    << arg2 << ".");

  DeclException4(WrongNumberOfValues,
    std::string, std::string, unsigned int, unsigned int,
    << "In material " << arg1 << ", the number of values under the label "
    << arg2 << " (" << arg3 << ")"
    << " is inconsistent with number_of_groups = " << arg4 << ".");

  DeclException4(WrongSign,
    std::string, std::string, double, unsigned int,
    << "Material " << arg1 << " contains a value of the wrong sign in "
    << arg2 << "." << " (" << arg3 << " at index " << arg4 << ")");

  DeclException4(EnergyGroupBoundariesNotSorted,
    std::string, double, double, double,
    << "Energy group boundary vector for material " << arg1
    << " contains {" << arg2 << ", " << arg3 << ", "
    << arg4 << "}, which are out of order.");

  /*!
    CheckValid expects normalization to machine precision levels,
    which means experimental Chi data with lower precision would
    have to be changed before going into the .material file
  */
  DeclException2(ChiDoesNotSumToOne, std::string, double,
    << "In material " << arg1 << ", the sum of values in the Chi vector is "
    << std::setprecision(std::numeric_limits<double>::digits10 + 1)
    << arg2 << ". It must be normalized to 1.");
};

#endif // BART_SRC_MATERIAL_MATERIAL_PROPERTIES_H_
