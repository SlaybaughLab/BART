#include "data/system/right_hand_side.h"

#include <utility>

namespace bart {

namespace data {

namespace system {

RightHandSide::RightHandSide(std::unordered_set<VariableTerms> variable_terms)
    : variable_terms_(variable_terms)
{
  for (auto term : variable_terms) {
    RightHandSidePtrMap variable_term_map;
    variable_right_hand_side_terms_[term] = variable_term_map;
  }
}


void RightHandSide::SetFixedTermPtr(Index index, std::shared_ptr<MPIVector> to_set) {
  fixed_right_hand_side_ptrs_[index] = to_set;
}

void RightHandSide::SetFixedTermPtr(GroupNumber group, std::shared_ptr<MPIVector> to_set) {
  SetFixedTermPtr({group, 0}, to_set);
}

std::shared_ptr<MPIVector> RightHandSide::GetFixedTermPtr(Index index) {
  try {
    return fixed_right_hand_side_ptrs_.at(index);
  } catch (std::out_of_range &exc) {
    return nullptr;
  }
}

std::shared_ptr<MPIVector> RightHandSide::GetFixedTermPtr(GroupNumber group) {
  return GetFixedTermPtr({group, 0});
}

void RightHandSide::SetVariableTermPtr(Index index,
                                       VariableTerms term,
                                       std::shared_ptr<MPIVector> to_set) {
  AssertThrow(variable_terms_.count(term) != 0,
              dealii::ExcMessage("Tried to set a right hand side with a variable "
                                 "term that it does not have set as variable"));
  variable_right_hand_side_terms_[term][index] = to_set;
}

void RightHandSide::SetVariableTermPtr(GroupNumber group,
                                       VariableTerms term,
                                       std::shared_ptr<MPIVector> to_set) {
  SetVariableTermPtr({group, 0}, term, to_set);
}

std::shared_ptr<MPIVector> RightHandSide::GetVariableTermPtr(Index index,
                                                             VariableTerms term) {
  AssertThrow(variable_terms_.count(term) != 0,
              dealii::ExcMessage("Tried to access a right hand side with a variable "
                                 "term that it does not have set as variable"));
  try {
    return variable_right_hand_side_terms_[term].at(index);
  } catch (std::out_of_range &exc) {
    return nullptr;
  }
}

std::shared_ptr<MPIVector> RightHandSide::GetVariableTermPtr(GroupNumber group,
                                                             VariableTerms term) {
  return GetVariableTermPtr({group, 0}, term);
}


} // namespace system

} // namespace data

} // namespace bart