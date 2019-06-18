#ifndef BART_SRC_DATA_SYSTEM_TERM_H_
#define BART_SRC_DATA_SYSTEM_TERM_H_

#include <memory>
#include <map>

#include "system/system_types.h"
#include "system/terms/term_i.h"

namespace bart {

namespace system {

namespace terms {
/*! \brief Stores and provides linear and bilinear terms.
 *
 * This class enables storage of any data storage type, but in general uses
 * matrices for bilinear terms, and vectors for linear terms.
 *
 * Data objects are stored differently if they are fixed (will not change iteration
 * to iteration) or variable (do change iteration to iteration). A single one
 * is stored for the fixed terms, as they can all be stamped once and left. Each
 * variable term has its own data object, so that they can be individually updated
 * when needed.
 *
 * Each group and angle requires its own data object, so storage is based on
 * an index, comprised of the group number, and the _index_ of the angle.
 * Overloads are provided that only require group number, which will retrieve
 * and store those with an angle index of zero.
 *
 * @tparam TermPair a std::pair that includes two types, (1) the storage type
 * and (2) an enum holding the possible variable terms.
 *
 * Example: If we wanted to create an MPI linear term (using vectors) that can
 * hold a fixed term and a variable term called kNewVariableTerm, we would use
 * the class as follows:
 *
 * \code{.cpp}
 *
 * using system::MPIVector;
 * enum class VariableTerms = { kNewVariableTerm };
 *
 * system::terms::Term<MPIVector, VariableTerms> new_term({kNewVariableTerm});
 *
 * \endcode
 *
 * You would then be able to set fixed vectors, and variable vectors for the
 * kNewVariableTerm term.
 *
 */
template <typename TermPair>
class Term : public TermI<TermPair> {
 public:
  using StorageType = typename TermPair::first_type;
  using VariableTermType = typename TermPair::second_type;

  /*! Constructor.
   *
   * When constructed, a set that indicates which members of VaraibleTermType
   * (the second type passed in the TermPair template parameter) will vary from
   * iteration to iteration. If not specified, no extra storage objects will
   * be created for the term, and trying to access them will result in an error.
   *
   * @param [variable_terms] unordered set that indicates the members of
   *                         VariableTermType that will vary.
   */
  explicit Term(std::unordered_set<VariableTermType> variable_terms = {});
  virtual ~Term() = default;

  std::unordered_set<VariableTermType> GetVariableTerms() const override {
    return variable_terms_;
  };

  void SetFixedTermPtr(Index index, std::shared_ptr<StorageType> to_set) override;
  void SetFixedTermPtr(GroupNumber group, std::shared_ptr<StorageType> to_set) override;
  std::shared_ptr<StorageType> GetFixedTermPtr(Index index) override;
  std::shared_ptr<StorageType> GetFixedTermPtr(GroupNumber group) override;

  void SetVariableTermPtr(Index index,
                          VariableTermType term,
                          std::shared_ptr<StorageType> to_set) override;
  void SetVariableTermPtr(GroupNumber group,
                          VariableTermType term,
                          std::shared_ptr<StorageType> to_set) override;
  std::shared_ptr<StorageType> GetVariableTermPtr(Index index,
                                                  VariableTermType term) override;
  std::shared_ptr<StorageType> GetVariableTermPtr(GroupNumber group,
                                                  VariableTermType term) override;



 private:
  const std::unordered_set<VariableTermType> variable_terms_;

  using TermPtrMap = std::map<Index, std::shared_ptr<StorageType>>;

  TermPtrMap fixed_term_ptrs_;

  std::map<VariableTermType, TermPtrMap> variable_term_ptrs_;
};

using MPILinearTerm = Term<system::terms::MPILinearTermPair>;
using MPIBilinearTerm = Term<system::terms::MPIBilinearTermPair>;

} // namespace terms

} // namespace system

} // namespace bart

#endif // BART_SRC_DATA_SYSTEM_TERM_H_