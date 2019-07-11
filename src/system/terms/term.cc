#include "term.h"

#include <utility>

namespace bart {

namespace system {

namespace terms {

template <typename TermPair>
Term<TermPair>::Term(std::unordered_set<VariableTermType> variable_terms)
    : variable_terms_(variable_terms)
{
  for (auto term : variable_terms) {
    TermPtrMap variable_term_map;
    variable_term_ptrs_[term] = variable_term_map;
  }
}
template <typename TermPair>
void Term<TermPair>::SetFixedTermPtr(Index index, std::shared_ptr<StorageType> to_set) {
  fixed_term_ptrs_[index] = to_set;
}

template <typename TermPair>
void Term<TermPair>::SetFixedTermPtr(GroupNumber group, std::shared_ptr<StorageType> to_set) {
  SetFixedTermPtr({group, 0}, to_set);
}

template <typename TermPair>
auto Term<TermPair>::GetFixedTermPtr(Index index) -> std::shared_ptr<StorageType> {
  try {
    return fixed_term_ptrs_.at(index);
  } catch (std::out_of_range &exc) {
    return nullptr;
  }
}

template <typename TermPair>
auto Term<TermPair>::GetFixedTermPtr(GroupNumber group) -> std::shared_ptr<StorageType> {
  return GetFixedTermPtr({group, 0});
}

template <typename TermPair>
void Term<TermPair>::SetVariableTermPtr(Index index,
                                       VariableTermType term,
                                       std::shared_ptr<StorageType> to_set) {
  AssertThrow(variable_terms_.count(term) != 0,
              dealii::ExcMessage("Tried to set a right hand side with a variable "
                                 "term that it does not have set as variable"));
  variable_term_ptrs_[term][index] = to_set;
}

template <typename TermPair>
void Term<TermPair>::SetVariableTermPtr(GroupNumber group,
                                       VariableTermType term,
                                       std::shared_ptr<StorageType> to_set) {
  SetVariableTermPtr({group, 0}, term, to_set);
}

template <typename TermPair>
auto Term<TermPair>::GetVariableTermPtr(Index index,
                                        VariableTermType term) -> std::shared_ptr<StorageType> {
  AssertThrow(variable_terms_.count(term) != 0,
              dealii::ExcMessage("Tried to access a right hand side with a variable "
                                 "term that it does not have set as variable"));
  try {
    return variable_term_ptrs_[term].at(index);
  } catch (std::out_of_range &exc) {
    return nullptr;
  }
}

template <typename TermPair>
auto Term<TermPair>::GetVariableTermPtr(GroupNumber group,
                                        VariableTermType term) -> std::shared_ptr<StorageType> {
  return GetVariableTermPtr({group, 0}, term);
}

template <>
std::shared_ptr<system::MPIVector> Term<MPILinearTermPair>::GetFullTermPtr(
    Index index) const {
  return std::shared_ptr<system::MPIVector>();
}

template <>
std::shared_ptr<system::MPISparseMatrix> Term<MPIBilinearTermPair>::GetFullTermPtr(
    Index index) const {
  return std::shared_ptr<system::MPISparseMatrix>();
}

template class Term<system::terms::MPILinearTermPair>;
template class Term<system::terms::MPIBilinearTermPair>;

} // namespace terms

} // namespace system

} // namespace bart