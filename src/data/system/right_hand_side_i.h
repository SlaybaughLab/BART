#ifndef BART_SRC_DATA_SYSTEM_RIGHT_HAND_SIDE_I_H_
#define BART_SRC_DATA_SYSTEM_RIGHT_HAND_SIDE_I_H_

#include <unordered_set>

#include "data/system/rhs_lhs_types.h"

namespace bart {

namespace data {

namespace system {

class RightHandSideI {
 public:
  enum class VariableTerms{
    kScatteringSource = 0,
    kFissionSource = 1,
  };
  virtual ~RightHandSideI() = default;


  virtual std::unordered_set<VariableTerms> GetVariableTerms() const = 0;

  virtual void SetFixedPtr(Index index, std::shared_ptr<MPIVector> to_set) = 0;
  virtual void SetFixedPtr(GroupNumber group, std::shared_ptr<MPIVector> to_set) = 0;
  virtual std::shared_ptr<MPIVector> GetFixedPtr(Index index) = 0;
  virtual std::shared_ptr<MPIVector> GetFixedPtr(GroupNumber group) = 0;

  virtual void SetVariablePtr(Index index,
                              VariableTerms term,
                              std::shared_ptr<MPIVector> to_set) = 0;
  virtual void SetVariablePtr(GroupNumber group,
                              VariableTerms term,
                              std::shared_ptr<MPIVector> to_set) = 0;
  virtual std::shared_ptr<MPIVector> GetVariablePtr(Index index,
                                                    VariableTerms term) = 0;

};

} // namespace system

} // namespace data

} // namespace bart

#endif // BART_SRC_DATA_SYSTEM_RIGHT_HAND_SIDE_I_H_