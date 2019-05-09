#ifndef BART_SRC_DATA_SYSTEM_RIGHT_HAND_SIDE_FIXED_H_
#define BART_SRC_DATA_SYSTEM_RIGHT_HAND_SIDE_FIXED_H_

#include <memory>
#include <map>

#include "data/system/system_types.h"
#include "data/system/right_hand_side_i.h"

namespace bart {

namespace data {

namespace system {

class RightHandSide : public RightHandSideI {
 public:
  explicit RightHandSide(std::unordered_set<VariableTerms> = {});
  virtual ~RightHandSide() = default;

  std::unordered_set<VariableTerms> GetVariableTerms() const override {
    return variable_terms_;
  };

  void SetFixedTermPtr(Index index, std::shared_ptr<MPIVector> to_set) override;
  void SetFixedTermPtr(GroupNumber group, std::shared_ptr<MPIVector> to_set) override;
  std::shared_ptr<MPIVector> GetFixedTermPtr(Index index) override;
  std::shared_ptr<MPIVector> GetFixedTermPtr(GroupNumber group) override;

  void SetVariableTermPtr(Index index,
                          VariableTerms term,
                          std::shared_ptr<MPIVector> to_set) override;
  void SetVariableTermPtr(GroupNumber group,
                          VariableTerms term,
                          std::shared_ptr<MPIVector> to_set) override;
  std::shared_ptr<MPIVector> GetVariableTermPtr(Index index,
                                                VariableTerms term) override;
  std::shared_ptr<MPIVector> GetVariableTermPtr(GroupNumber group,
                                                VariableTerms term) override;



 private:
  const std::unordered_set<VariableTerms> variable_terms_;

  using RightHandSidePtrMap = std::map<Index, std::shared_ptr<MPIVector>>;

  RightHandSidePtrMap fixed_right_hand_side_ptrs_;

  std::map<VariableTerms, RightHandSidePtrMap> variable_right_hand_side_terms_;
};


} // namespace system

} // namespace data

} // namespace bart

#endif // BART_SRC_DATA_SYSTEM_RIGHT_HAND_SIDE_FIXED_H_