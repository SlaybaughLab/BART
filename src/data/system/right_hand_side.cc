#include "data/system/right_hand_side.h"

namespace bart {

namespace data {

namespace system {

RightHandSide::RightHandSide(std::unordered_set<VariableTerms> variable_terms)
    : variable_terms_(variable_terms)
{}



} // namespace system

} // namespace data

} // namespace bart