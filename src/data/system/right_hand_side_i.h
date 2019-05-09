#ifndef BART_SRC_DATA_SYSTEM_RIGHT_HAND_SIDE_I_H_
#define BART_SRC_DATA_SYSTEM_RIGHT_HAND_SIDE_I_H_

#include <unordered_set>

#include "data/system/system_types.h"

namespace bart {

namespace data {

namespace system {

/*! Stores and provides right hand side vectors.
 *
 * Vectors are stored differently if they are fixed (will not change iteration
 * to iteration) or variable (do change iteration to iteration). A single vector
 * is stored for the fixed terms, as they can all be stamped once and left. Each
 * variable term has its own vector, so that they can be individually updated
 * when needed.
 *
 * Each group and angle requires its own vector, so vectors are stored based on
 * an index, comprised of the group number, and the _index_ of the angle.
 * Overloads are provided that only require group number, which will retrieve
 * and store vectors with an angle index of zero.
 */
class RightHandSideI {
 public:
  //! Variable terms that can be stored
  enum class VariableTerms{
    kScatteringSource = 0, //!< Scattering source
    kFissionSource = 1,    //!< Fission source
  };
  virtual ~RightHandSideI() = default;

  /*! \brief Returns the set of terms that are set as variable. */
  virtual std::unordered_set<VariableTerms> GetVariableTerms() const = 0;

  /*! \brief Sets the MPIVector (a pointer to the vector) used to store the fixed term.
   *
   * @param index corresponding index for the MPIVector
   * @param to_set pointer to the vector to set
   */
  virtual void SetFixedTermPtr(Index index, std::shared_ptr<MPIVector> to_set) = 0;
  /*! \brief Sets the MPIVector (a pointer to the vector) used to store the fixed term.
   *
   * @param group group number to be set
   * @param to_set pointer to the vector to set
   */
  virtual void SetFixedTermPtr(GroupNumber group, std::shared_ptr<MPIVector> to_set) = 0;
  /*! \brief Returns a pointer to the fixed term for a given index
   *
   * @param index index for the fixed term
   * @return a pointer to the fixed term
   */
  virtual std::shared_ptr<MPIVector> GetFixedTermPtr(Index index) = 0;
  /*! \brief Returns a pointer to the fixed term for a given group.
   *
   * @param group group to return
   * @return a pointer to the fixed term
   */
  virtual std::shared_ptr<MPIVector> GetFixedTermPtr(GroupNumber group) = 0;

  /*! \brief Sets the MPIVector (a pointer to the vector) used to store the
   *         specified fixed term.
   *
   * @param index index to be set
   * @param term specifies the variable term to be set
   * @param to_set pointer to the vector to be set
   */
  virtual void SetVariableTermPtr(Index index,
                                  VariableTerms term,
                                  std::shared_ptr<MPIVector> to_set) = 0;
  /*! \brief Sets the MPIVector (a pointer to the vector) used to store the
   *         specified fixed term.
   *
   * @param group group number to be set
   * @param term specifies the variable term to be set
   * @param to_set pointer to the vector to be set
   */
  virtual void SetVariableTermPtr(GroupNumber group,
                                  VariableTerms term,
                                  std::shared_ptr<MPIVector> to_set) = 0;
  /*! \brief Returns a pointer to the specified variable term.
   *
   * @param index index for the variable term
   * @param term variable term to return
   * @return pointer to the variable term.
   */
  virtual std::shared_ptr<MPIVector> GetVariableTermPtr(Index index,
                                                        VariableTerms term) = 0;
  /*! \brief Returns a pointer to the specified variable term.
   *
   * @param group group to return
   * @param term variable term to return
   * @return pointer to the variable term.
   */
  virtual std::shared_ptr<MPIVector> GetVariableTermPtr(GroupNumber group,
                                                        VariableTerms term) = 0;

};

} // namespace system

} // namespace data

} // namespace bart

#endif // BART_SRC_DATA_SYSTEM_RIGHT_HAND_SIDE_I_H_