#ifndef BART_SRC_DATA_SYSTEM_TESTS_BILINEAR_TERM_MOCK_HPP_
#define BART_SRC_DATA_SYSTEM_TESTS_BILINEAR_TERM_MOCK_HPP_

#include "system/terms/term_i.h"

#include <unordered_set>

#include "system/system_types.h"
#include "system/terms/term_types.h"
#include "test_helpers/gmock_wrapper.h"

namespace bart::system::terms {

class BilinearTermMock : public TermI<system::terms::MPIBilinearTermPair> {
 public:
  using MPISparseMatrix = bart::system::MPISparseMatrix;
  using VariableBilinearTerms = system::terms::VariableBilinearTerms;
  MOCK_METHOD(std::shared_ptr<MPISparseMatrix>, GetFixedTermPtr, (Index), (override));
  MOCK_METHOD(std::shared_ptr<MPISparseMatrix>, GetFixedTermPtr, (GroupNumber), (override));
  MOCK_METHOD(std::shared_ptr<MPISparseMatrix>, GetFullTermPtr, (Index), (const, override));
  MOCK_METHOD(std::shared_ptr<MPISparseMatrix>, GetVariableTermPtr, (Index, VariableBilinearTerms), (override));
  MOCK_METHOD(std::shared_ptr<MPISparseMatrix>, GetVariableTermPtr, (GroupNumber, VariableBilinearTerms), (override));
  MOCK_METHOD(std::unordered_set<VariableBilinearTerms>, GetVariableTerms, (), (const, override));
  MOCK_METHOD(void, SetFixedTermPtr, (Index, std::shared_ptr<MPISparseMatrix>), (override));
  MOCK_METHOD(void, SetFixedTermPtr, (GroupNumber, std::shared_ptr<MPISparseMatrix>), (override));
  MOCK_METHOD(void, SetVariableTermPtr, (Index, VariableBilinearTerms, std::shared_ptr<MPISparseMatrix>), (override));
  MOCK_METHOD(void, SetVariableTermPtr, (GroupNumber, VariableBilinearTerms, std::shared_ptr<MPISparseMatrix>), (override));
};

} // namespace bart::system::terms

#endif // BART_SRC_DATA_SYSTEM_TESTS_BILINEAR_TERM_MOCK_HPP_