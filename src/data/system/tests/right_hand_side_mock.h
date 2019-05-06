#ifndef BART_SRC_DATA_SYSTEM_TESTS_RIGHT_HAND_SIDE_MOCK_H_
#define BART_SRC_DATA_SYSTEM_TESTS_RIGHT_HAND_SIDE_MOCK_H_

#include <unordered_set>

#include "data/system/right_hand_side_i.h"
#include "data/system/system_types.h"

namespace bart {

namespace data {

namespace system {

class RightHandSideMock : public RightHandSideI {
 public:
  using RightHandSideI::VariableTerms;

  MOCK_CONST_METHOD0(GetVariableTerms, std::unordered_set<VariableTerms>());
  MOCK_METHOD2(SetFixedPtr, void(Index, std::shared_ptr<MPIVector>));
  MOCK_METHOD2(SetFixedPtr, void(GroupNumber, std::shared_ptr<MPIVector>));
  MOCK_METHOD1(GetFixedPtr, std::shared_ptr<MPIVector>(Index));
  MOCK_METHOD1(GetFixedPtr, std::shared_ptr<MPIVector>(GroupNumber));
  MOCK_METHOD3(SetVariablePtr, void(Index, VariableTerms,
      std::shared_ptr<MPIVector>));
  MOCK_METHOD3(SetVariablePtr, void(GroupNumber, VariableTerms,
      std::shared_ptr<MPIVector>));
  MOCK_METHOD2(GetVariablePtr, std::shared_ptr<MPIVector>(Index, VariableTerms));
  MOCK_METHOD2(GetVariablePtr, std::shared_ptr<MPIVector>(GroupNumber, VariableTerms));
};

} // namespace system

} // namespace data

} // namespace bart

#endif // BART_SRC_DATA_SYSTEM_TESTS_RIGHT_HAND_SIDE_MOCK_H_