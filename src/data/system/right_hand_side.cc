#include "data/system/right_hand_side.h"

namespace bart {

namespace data {

namespace system {

RightHandSide::RightHandSide(std::unordered_set<VariableTerms> variable_terms)
    : variable_terms_(variable_terms)
{}


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

}

} // namespace system

} // namespace data

} // namespace bart