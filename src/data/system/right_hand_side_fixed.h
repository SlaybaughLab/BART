#ifndef BART_SRC_DATA_SYSTEM_RIGHT_HAND_SIDE_FIXED_H_
#define BART_SRC_DATA_SYSTEM_RIGHT_HAND_SIDE_FIXED_H_

#include <memory>

#include "data/system/rhs_lhs_types.h"
#include "data/system/right_hand_side_i.h"

namespace bart {

namespace data {

namespace system {

class RightHandSideFixed : public RightHandSideI {
 public:
  virtual ~RightHandSideFixed() = default;

  virtual std::shared_ptr<MPIVector> GetFixedPtr(Index index) {
    try {
      return fixed_right_hand_side_.at(index);
    } catch (std::out_of_range &exc) {
      return nullptr;
    }
  };

  virtual std::shared_ptr<MPIVector> GetFixedPtr(GroupNumber group) {
    return GetFixedPtr({group, 0});
  };

  void SetFixedPtr(Index index, std::shared_ptr<MPIVector> to_set) {
    fixed_right_hand_side_[index] = to_set;
  };

  void SetFixedPtr(GroupNumber group, std::shared_ptr<MPIVector> to_set) {
    SetFixedPtr({group, 0}, to_set);
  }


 private:
  std::map<Index, std::shared_ptr<MPIVector>> fixed_right_hand_side_;
};


} // namespace system

} // namespace data

} // namespace bart

#endif // BART_SRC_DATA_SYSTEM_RIGHT_HAND_SIDE_FIXED_H_