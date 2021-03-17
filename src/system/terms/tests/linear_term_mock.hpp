#ifndef BART_SRC_DATA_SYSTEM_TESTS_LINEAR_TERM_MOCK_HPP_
#define BART_SRC_DATA_SYSTEM_TESTS_LINEAR_TERM_MOCK_HPP_

#include <unordered_set>

#include "system/system_types.h"
#include "system/terms/term_i.h"
#include "system/terms/term_types.h"

namespace bart::system::terms {

class LinearTermMock : public TermI<system::terms::MPILinearTermPair> {
 public:
  using MPIVector = bart::system::MPIVector;
  using VariableLinearTerms = system::terms::VariableLinearTerms;
  MOCK_METHOD(std::shared_ptr<MPIVector>, GetFixedTermPtr, (Index), (override));
  MOCK_METHOD(std::shared_ptr<MPIVector>, GetFixedTermPtr, (GroupNumber));
  MOCK_METHOD(std::shared_ptr<MPIVector>, GetFullTermPtr, (Index), (const, override));
  MOCK_METHOD(std::shared_ptr<MPIVector>, GetVariableTermPtr, (Index, VariableLinearTerms), (override));
  MOCK_METHOD(std::shared_ptr<MPIVector>, GetVariableTermPtr, (GroupNumber, VariableLinearTerms), (override));
  MOCK_METHOD(std::unordered_set<VariableLinearTerms>, GetVariableTerms, (), (const, override));
  MOCK_METHOD(void, SetFixedTermPtr, (Index, std::shared_ptr<MPIVector>), (override));
  MOCK_METHOD(void, SetFixedTermPtr, (GroupNumber, std::shared_ptr<MPIVector>), (override));
  MOCK_METHOD(void, SetVariableTermPtr, (Index, VariableLinearTerms, std::shared_ptr<MPIVector>), (override));
  MOCK_METHOD(void, SetVariableTermPtr, (GroupNumber, VariableLinearTerms, std::shared_ptr<MPIVector>), (override));
};

} // namespace bart::system::terms

#endif // BART_SRC_DATA_SYSTEM_TESTS_LINEAR_TERM_MOCK_HPP_