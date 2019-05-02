#ifndef BART_SRC_DATA_SYSTEM_RIGHT_HAND_SIDE_I_H_
#define BART_SRC_DATA_SYSTEM_RIGHT_HAND_SIDE_I_H_

#include "data/system/rhs_lhs_types.h"

namespace bart {

namespace data {

namespace system {

class RightHandSideI {
 public:
  virtual ~RightHandSideI() = default;
  virtual std::shared_ptr<MPIVector> GetFixedPtr(Index index) = 0;
  virtual std::shared_ptr<MPIVector> GetFixedPtr(GroupNumber group) = 0;
  virtual void SetFixedPtr(Index index, std::shared_ptr<MPIVector> to_set) = 0;
  virtual void SetFixedPtr(GroupNumber group, std::shared_ptr<MPIVector> to_set) = 0;
};

} // namespace system

} // namespace data

} // namespace bart

#endif // BART_SRC_DATA_SYSTEM_RIGHT_HAND_SIDE_I_H_