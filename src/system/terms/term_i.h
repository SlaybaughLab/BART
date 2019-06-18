#ifndef BART_SRC_DATA_SYSTEM_TERM_I_H_
#define BART_SRC_DATA_SYSTEM_TERM_I_H_

#include <unordered_set>

#include "system/terms/term_types.h"
#include "system/system_types.h"

namespace bart {

namespace system {

namespace terms {

template <typename TermPair>
class TermI {
 public:
  using StorageType = typename TermPair::first_type;
  using VariableTermType = typename TermPair::second_type;
  virtual ~TermI() = default;

  /*! \brief Returns the set of terms that are set as variable. */
  virtual std::unordered_set<VariableTermType> GetVariableTerms() const = 0;

  /*! \brief Sets a pointer to the data object used to store the fixed term.
   *
   * @param index corresponding index for the data object
   * @param to_set pointer to the data object to set
   */
  virtual void SetFixedTermPtr(Index index, std::shared_ptr<StorageType> to_set) = 0;
  /*! \brief Sets a pointer to the data object used to store the fixed term.
   *
   * @param group corresponding group number for the data object
   * @param to_set pointer to the data object to set
   */
  virtual void SetFixedTermPtr(GroupNumber group, std::shared_ptr<StorageType> to_set) = 0;
  /*! \brief Returns a pointer to the fixed term for a given index
   *
   * @param index index for the fixed term
   * @return a pointer to the fixed term
   */
  virtual std::shared_ptr<StorageType> GetFixedTermPtr(Index index) = 0;
  /*! \brief Returns a pointer to the fixed term for a given group.
   *
   * @param group group to return
   * @return a pointer to the fixed term
   */
  virtual std::shared_ptr<StorageType> GetFixedTermPtr(GroupNumber group) = 0;

  /*! \brief Sets the MPIVector (a pointer to the data object) used to store the
   *         specified fixed term.
   *
   * @param index index to be set
   * @param term specifies the variable term to be set
   * @param to_set pointer to the data object to be set
   */
  virtual void SetVariableTermPtr(Index index,
                                  VariableTermType term,
                                  std::shared_ptr<StorageType> to_set) = 0;
  /*! \brief Sets the MPIVector (a pointer to the data object) used to store the
   *         specified fixed term.
   *
   * @param group group number to be set
   * @param term specifies the variable term to be set
   * @param to_set pointer to the data object to be set
   */
  virtual void SetVariableTermPtr(GroupNumber group,
                                  VariableTermType term,
                                  std::shared_ptr<StorageType> to_set) = 0;
  /*! \brief Returns a pointer to the specified variable term.
   *
   * @param index index for the variable term
   * @param term variable term to return
   * @return pointer to the variable term.
   */
  virtual std::shared_ptr<StorageType> GetVariableTermPtr(Index index,
                                                          VariableTermType term) = 0;
  /*! \brief Returns a pointer to the specified variable term.
   *
   * @param group group to return
   * @param term variable term to return
   * @return pointer to the variable term.
   */
  virtual std::shared_ptr<StorageType> GetVariableTermPtr(GroupNumber group,
                                                          VariableTermType term) = 0;

};

using MPILinearTermI = TermI<system::terms::MPILinearTermPair>;
using MPIBilinearTermI = TermI<system::terms::MPIBilinearTermPair>;

} // namespace terms

} // namespace system

} // namespace bart

#endif // BART_SRC_DATA_SYSTEM_TERM_I_H_