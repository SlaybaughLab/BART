#ifndef BART_SRC_DATA_SYSTEM_TERM_H_
#define BART_SRC_DATA_SYSTEM_TERM_H_

#include <memory>
#include <map>

#include "data/system/system_types.h"
#include "data/system/term_i.h"

namespace bart {

namespace data {

namespace system {

template <typename TermPair>
class Term : public TermI<TermPair> {
 public:
  using StorageType = typename TermPair::first_type;
  using VariableTermType = typename TermPair::second_type;
//  explicit Term(std::unordered_set<VariableTermType> = {});
//  virtual ~Term() = default;
//
//  std::unordered_set<VariableTermType> GetVariableTerms() const override {
//    return variable_terms_;
//  };
//
//  void SetFixedTermPtr(Index index, std::shared_ptr<StorageType> to_set) override;
//  void SetFixedTermPtr(GroupNumber group, std::shared_ptr<StorageType> to_set) override;
//  std::shared_ptr<StorageType> GetFixedTermPtr(Index index) override;
//  std::shared_ptr<StorageType> GetFixedTermPtr(GroupNumber group) override;
//
//  void SetVariableTermPtr(Index index,
//                          VariableTermType term,
//                          std::shared_ptr<StorageType> to_set) override;
//  void SetVariableTermPtr(GroupNumber group,
//                          VariableTermType term,
//                          std::shared_ptr<StorageType> to_set) override;
//  std::shared_ptr<StorageType> GetVariableTermPtr(Index index,
//                                                VariableTermType term) override;
//  std::shared_ptr<StorageType> GetVariableTermPtr(GroupNumber group,
//                                                VariableTermType term) override;
//
//
//
// private:
//  const std::unordered_set<VariableTermType> variable_terms_;
//
//  using TermPtrMap = std::map<Index, std::shared_ptr<StorageType>>;
//
//  TermPtrMap fixed_right_hand_side_ptrs_;
//
//  std::map<VariableTermType, TermPtrMap> variable_right_hand_side_terms_;
};


} // namespace system

} // namespace data

} // namespace bart

#endif // BART_SRC_DATA_SYSTEM_TERM_H_