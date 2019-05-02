#ifndef BART_SRC_DATA_SYSTEM_RIGHT_HAND_SIDE_FIXED_H_
#define BART_SRC_DATA_SYSTEM_RIGHT_HAND_SIDE_FIXED_H_

#include <memory>
#include <vector>

#include "data/system/rhs_lhs_types.h"
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

  void SetFixedPtr(Index index, std::shared_ptr<MPIVector> to_set) override;
  void SetFixedPtr(GroupNumber group, std::shared_ptr<MPIVector> to_set) override;
  std::shared_ptr<MPIVector> GetFixedPtr(Index index) override;
  std::shared_ptr<MPIVector> GetFixedPtr(GroupNumber group) override;

  void SetVariablePtr(Index index,
                      VariableTerms term,
                      std::shared_ptr<MPIVector> to_set) override;

 private:
  const std::unordered_set<VariableTerms> variable_terms_;
  std::map<Index, std::shared_ptr<MPIVector>> fixed_right_hand_side_ptrs_;
};


} // namespace system

} // namespace data

} // namespace bart

#endif // BART_SRC_DATA_SYSTEM_RIGHT_HAND_SIDE_FIXED_H_