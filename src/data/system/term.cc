#include "data/system/term.h"

#include <utility>

namespace bart {

namespace data {

namespace system {

//template <typename StorageType, typename VariableTermType>
//Term::Term(std::unordered_set<VariableTermType> variable_terms)
//    : variable_terms_(variable_terms)
//{
//  for (auto term : variable_terms) {
//    TermPtrMap variable_term_map;
//    variable_right_hand_side_terms_[term] = variable_term_map;
//  }
//}
//
//template <typename StorageType, typename VariableTermType>
//void Term::SetFixedTermPtr(Index index, std::shared_ptr<StorageType> to_set) {
//  fixed_right_hand_side_ptrs_[index] = to_set;
//}
//
//template <typename StorageType, typename VariableTermType>
//void Term::SetFixedTermPtr(GroupNumber group, std::shared_ptr<StorageType> to_set) {
//  SetFixedTermPtr({group, 0}, to_set);
//}
//
//template <typename StorageType, typename VariableTermType>
//std::shared_ptr<StorageType> Term::GetFixedTermPtr(Index index) {
//  try {
//    return fixed_right_hand_side_ptrs_.at(index);
//  } catch (std::out_of_range &exc) {
//    return nullptr;
//  }
//}
//
//template <typename StorageType, typename VariableTermType>
//std::shared_ptr<StorageType> Term::GetFixedTermPtr(GroupNumber group) {
//  return GetFixedTermPtr({group, 0});
//}
//
//template <typename StorageType, typename VariableTermType>
//void Term::SetVariableTermPtr(Index index,
//                                       VariableTermType term,
//                                       std::shared_ptr<StorageType> to_set) {
//  AssertThrow(variable_terms_.count(term) != 0,
//              dealii::ExcMessage("Tried to set a right hand side with a variable "
//                                 "term that it does not have set as variable"));
//  variable_right_hand_side_terms_[term][index] = to_set;
//}
//
//template <typename StorageType, typename VariableTermType>
//void Term::SetVariableTermPtr(GroupNumber group,
//                                       VariableTermType term,
//                                       std::shared_ptr<StorageType> to_set) {
//  SetVariableTermPtr({group, 0}, term, to_set);
//}
//
//template <typename StorageType, typename VariableTermType>
//std::shared_ptr<StorageType> Term::GetVariableTermPtr(Index index,
//                                                             VariableTermType term) {
//  AssertThrow(variable_terms_.count(term) != 0,
//              dealii::ExcMessage("Tried to access a right hand side with a variable "
//                                 "term that it does not have set as variable"));
//  try {
//    return variable_right_hand_side_terms_[term].at(index);
//  } catch (std::out_of_range &exc) {
//    return nullptr;
//  }
//}
//
//template <typename StorageType, typename VariableTermType>
//std::shared_ptr<StorageType> Term::GetVariableTermPtr(GroupNumber group,
//                                                             VariableTermType term) {
//  return GetVariableTermPtr({group, 0}, term);
//}


} // namespace system

} // namespace data

} // namespace bart