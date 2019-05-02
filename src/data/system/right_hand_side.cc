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


void RightHandSide::SetFixedPtr(Index index, std::shared_ptr<MPIVector> to_set) {
  fixed_right_hand_side_ptrs_[index] = to_set;
}

void RightHandSide::SetFixedPtr(GroupNumber group, std::shared_ptr<MPIVector> to_set) {
  SetFixedPtr({group, 0}, to_set);
}

std::shared_ptr<MPIVector> RightHandSide::GetFixedPtr(Index index) {
  try {
    return fixed_right_hand_side_ptrs_.at(index);
  } catch (std::out_of_range &exc) {
    return nullptr;
  }
}

std::shared_ptr<MPIVector> RightHandSide::GetFixedPtr(GroupNumber group) {
  return GetFixedPtr({group, 0});
}

void RightHandSide::SetVariablePtr(Index index,
                                   VariableTerms term,
                                   std::shared_ptr<MPIVector> to_set) {
  AssertThrow(variable_terms_.count(term) != 0,
              dealii::ExcMessage("Tried to set a right hand side with a variable "
                                 "term that it does not have set as variable"));
  variable_right_hand_side_terms_[term][index] = to_set;
}

void RightHandSide::SetVariablePtr(GroupNumber group,
                                   VariableTerms term,
                                   std::shared_ptr<MPIVector> to_set) {
  SetVariablePtr({group, 0}, term, to_set);
}



} // namespace system

} // namespace data

} // namespace bart