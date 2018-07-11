#ifndef BART_SRC_MATERIAL_MATERIAL_PROPERTIES_H_
#define BART_SRC_MATERIAL_MATERIAL_PROPERTIES_H_

#include <string>
#include <unordered_map>
#include <vector>
#include <unordered_set>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/lac/full_matrix.h>

#include "../../material.pb.h"

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
    constructor using strings from .material files
   */
  MaterialProperties(std::unordered_map<int, std::string> materials_as_strings /*! can be serialized or human-readable protobuf strings */,
                      bool is_eigen_problem,
                      bool do_nda /*! currently unused, may be removed in the future */,
                      int number_of_groups,
                      int number_of_materials,
                      std::unordered_set<int> fissile_ids = {} /*! only used in eigen problems */);

  /*!
    gets the necessary information from the parameter handler and delegates to the other constructor
  */
  MaterialProperties(dealii::ParameterHandler& prm);

  //! Class destructor.
  ~MaterialProperties();

  /*!
   A function to retrieve mapping: material id->if material is fissile boolean.

   \return A Hash table representing the desirable mapping.
   */
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
  std::unordered_map<int, std::vector<double>> GetNuSigf() const;

  //! A function to retrieve all scattering transfer matrices.
  std::unordered_map<int, dealii::FullMatrix<double>> GetSigS() const;

  //! A function to retrieve all scattering transfer matrices scaled by \f$4\pi\f$.
  std::unordered_map<int, dealii::FullMatrix<double>> GetSigSPerSter() const;

  //! A function to retrieve \f$\chi\nu\sigma_\mathrm{f}\f$ for all fissile materials.
  std::unordered_map<int, dealii::FullMatrix<double>> GetChiNuSigF() const;

  //! A function to retrieve \f$\chi\nu\sigma_\mathrm{f}/(4\pi)\f$ for all fissile materials.
  std::unordered_map<int, dealii::FullMatrix<double>> GetChiNuSigFPerSter() const;

  /*!
    tests if a Material is self-consistent and contains the required data

    the energy groups, sigma_t, and the scattering matrix are always required
    nu Sigma_f and Chi (normalized to 1) are required for fissile materials
  */
  static void CheckValid(const Material& material, bool require_fission_data = false);

  /*!
    tests if the given map holds materials that are all compatible with one another
    will maybe print a warning without failing if multiple IDs map to one material
      is printing a warning kosher?
    the following things are required to be the same in all of the materials:
      number of groups
      energy group boundaries
  */
  static void CheckConsistent(const std::unordered_map<int, Material>& materials);


 private:

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

  //! protobuf generated Material object for each material ID
  std::unordered_map<int, Material> materials_;

  //! Reads .material files listed in the parameter handler into strings
  static std::unordered_map<int, std::string> ReadMaterialStrings(dealii::ParameterHandler& prm);

  //! Reads fissile material ids listed in the parameter handler
  static std::unordered_set<int> ReadFissileIDs(dealii::ParameterHandler& prm);

  /*!
    populates std::unordered_map<int, Material> materials_ from strings
    accepts both serialized and human-readable protobuf files consistent with material.proto
  */
  void ParseMaterials(const std::unordered_map<int, std::string>& materials_as_strings);

  //! populates is_material_fissile_ with all material IDs and a bool if it is in fissile_ids_
  void PopulateFissileMap();

  //! populates all Get-able material data from materials_
  void PopulateData();

  //! throws an exception if fissile_ids_ contains an unused material id or if it is empty in an eigen problem
  void CheckFissileIDs() const;

  //! throws an exception if n_material_ doesn't equal the number of materials read in
  void CheckNumberOfMaterials() const;

  //! throws an exception if n_group_ doesn't equal the number of groups in each material
  void CheckNumberOfGroups() const;

  //! calls checkValid on each material
  void CheckValidEach() const;

  //! calls CheckConsistent on materials_;
  void CheckConsistent() const;

  //! returns full name, abbreviation, and id of a material into one string
  static std::string CombinedName(const Material& material, const std::string& delimiter = "|");

 public:
  /*
   * dealII macros for declaring exceptions are documented in https://www.dealii.org/9.0.0/doxygen/deal.II/group__Exceptions.html
   * should these be public?
   * using exceptions goes againt the Google style guide, but under these circumstances we want the program to terminate
   * check that this is kosher
  */

  DeclExceptionMsg(NoFissileIDs, "At least one material ID must be specified as fissile for eigen problems.");
  
  DeclException1(FissileIDInvalid, int,
    << "Material ID " << arg1 << " was specified as fissile, but no material with ID " << arg1 << " exists.");

  DeclException1(FailedToParseMaterialString, int,
    << "Could not parse protobuf input for material ID " << arg1 << " as serialized or human-readable format.");
  
  DeclException2(WrongNumberOfMaterials, int, unsigned int,
    << "The actual number of materials read in (" << arg2
    << ") does not match the number of materials specified (" << arg1 << ").");

  DeclException4(WrongNumberOfGroups, int, std::string, unsigned int, int,
    << "The number_of_groups in material number " << arg1
    << " (" << arg2 << ") does not match the number of groups in MaterialProperties."
    << " (" << arg3 << " != " << arg4 << ")");

  DeclException4(NumberOfGroupsMismatch, int, std::string, int, std::string,
    << "number of groups in material number " << arg1 << " (" << arg2 << ") and material number "
    << arg3 << " (" << arg4 << ") are different.");

  DeclException5(EnergyGroupsMismatch, int, std::string, int, std::string, unsigned int,
    << "number of groups in material number " << arg1 << " (" << arg2 << ") and material number "
    << arg3 << " (" << arg4 << ") differ at index " << arg5 << ".");

  // invalid material exceptions follow, TODO order these in the order that they're tested

  DeclException2(MissingProperty, std::string, std::string,
    << "Material with full_name|abbreviation|id  = " << arg1
    << " is missing required property " << arg2 << ".");

  DeclException2(MultipleDefinition, std::string, std::string,
    << "Material with full_name|abbreviation|id  = " << arg1
    << " has multiple definitions of " << arg2 << ", which is not allowed.");

  DeclException4(WrongNumberOfValues, std::string, std::string, unsigned int, unsigned int,
    << "In material with full_name|abbreviation|id  = " << arg1
    << ", the number of values under the label " << arg2 << " (" << arg3 << ")"
    << " is inconsistent with number_of_groups = " << arg4 << ".");

  /*!
   * CheckValid expects normalization to machine precision levels, which means experimental Chi data with
   * lower precision would have to be transformed before going into the .material file
   */
  DeclException2(ChiDoesNotSumToOne, std::string, double,
    << "In material with full_name|abbreviation|id  = " << arg1
    << ", the sum of values in the Chi vector is "
    << std::setprecision(std::numeric_limits<double>::digits10 + 1) << arg2 << ". It must be normalized to 1.");

  DeclException4(EnergyGroupBoundariesNotSorted, std::string, double, double, double,
    << "Energy group boundary vector for material with full_name|abbreviation|id  = " << arg1
    << " contains {" << arg2 << ", " << arg3 << ", " << arg4 << "}, which are out of order.");

  DeclException4(WrongSign, std::string, std::string, double, int,
    << "Material with full_name|abbreviation|id  = " << arg1
    << " contains value of the wrong sign in " << arg2 << "."
    << "(" << arg3 << " at index " << arg4 << ")");

  //! wrapper for for adding the material ID to an exception
  class ExceptionWithID : public dealii::ExceptionBase {
   public:
    ExceptionWithID(const dealii::ExceptionBase& exc, const int& id) : exc_(exc), id_(id) {};
    virtual ~ExceptionWithID() noexcept {}
    virtual void print_info(std::ostream &out) const {
      out << "Exception in material number " << id_ << ": ";
      exc_.print_info(out);
    }
    const dealii::ExceptionBase exc_;
    const int id_;
  };

};

#endif //BART_SRC_MATERIAL_MATERIAL_PROPERTIES_H_
