#ifndef BART_SRC_DATA_SYSTEM_TESTS_LINEAR_TERM_MOCK_H_
#define BART_SRC_DATA_SYSTEM_TESTS_LINEAR_TERM_MOCK_H_

#include <unordered_set>

#include "system/terms/term_i.h"
#include "system/system_types.h"

namespace bart {

namespace data {

namespace system {

class LinearTermMock : public TermI<data::system::MPILinearTermPair> {
 public:
  using MPIVector = bart::system::MPIVector;
  using VariableLinearTerms = data::system::VariableLinearTerms;

  MOCK_CONST_METHOD0(GetVariableTerms, std::unordered_set<VariableLinearTerms>());
  MOCK_METHOD2(SetFixedTermPtr, void(Index, std::shared_ptr<MPIVector>));
  MOCK_METHOD2(SetFixedTermPtr, void(GroupNumber, std::shared_ptr<MPIVector>));
  MOCK_METHOD1(GetFixedTermPtr, std::shared_ptr<MPIVector>(Index));
  MOCK_METHOD1(GetFixedTermPtr, std::shared_ptr<MPIVector>(GroupNumber));
  MOCK_METHOD3(SetVariableTermPtr, void(Index, VariableLinearTerms,
      std::shared_ptr<MPIVector>));
  MOCK_METHOD3(SetVariableTermPtr, void(GroupNumber, VariableLinearTerms,
      std::shared_ptr<MPIVector>));
  MOCK_METHOD2(GetVariableTermPtr, std::shared_ptr<MPIVector>(Index, VariableLinearTerms));
  MOCK_METHOD2(GetVariableTermPtr, std::shared_ptr<MPIVector>(GroupNumber, VariableLinearTerms));
};

} // namespace system

} // namespace data

} // namespace bart

#endif // BART_SRC_DATA_SYSTEM_TESTS_LINEAR_TERM_MOCK_H_