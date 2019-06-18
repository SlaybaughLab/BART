#ifndef BART_SRC_DATA_SYSTEM_TESTS_BILINEAR_TERM_MOCK_H_
#define BART_SRC_DATA_SYSTEM_TESTS_BILINEAR_TERM_MOCK_H_

#include "system/terms/term_i.h"

#include <unordered_set>

#include "system/system_types.h"
#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace data {

namespace system {

class BilinearTermMock : public TermI<data::system::MPIBilinearTermPair> {
 public:
  using MPISparseMatrix = bart::system::MPISparseMatrix;
  using VariableBilinearTerms = data::system::VariableBilinearTerms;

  MOCK_CONST_METHOD0(GetVariableTerms, std::unordered_set<VariableBilinearTerms>());
  MOCK_METHOD2(SetFixedTermPtr, void(Index, std::shared_ptr<MPISparseMatrix>));
  MOCK_METHOD2(SetFixedTermPtr, void(GroupNumber, std::shared_ptr<MPISparseMatrix>));
  MOCK_METHOD1(GetFixedTermPtr, std::shared_ptr<MPISparseMatrix>(Index));
  MOCK_METHOD1(GetFixedTermPtr, std::shared_ptr<MPISparseMatrix>(GroupNumber));
  MOCK_METHOD3(SetVariableTermPtr, void(Index, VariableBilinearTerms,
      std::shared_ptr<MPISparseMatrix>));
  MOCK_METHOD3(SetVariableTermPtr, void(GroupNumber, VariableBilinearTerms,
      std::shared_ptr<MPISparseMatrix>));
  MOCK_METHOD2(GetVariableTermPtr, std::shared_ptr<MPISparseMatrix>(Index, VariableBilinearTerms));
  MOCK_METHOD2(GetVariableTermPtr, std::shared_ptr<MPISparseMatrix>(GroupNumber, VariableBilinearTerms));
};

} // namespace system

} // namespace data

} // namespace bart

#endif // BART_SRC_DATA_SYSTEM_TESTS_BILINEAR_TERM_MOCK_H_