#ifndef BART_SRC_SYSTEM_SOLUTION_SOLUTION_H_
#define BART_SRC_SYSTEM_SOLUTION_SOLUTION_H_

#include "system/solution/solution_i.h"

namespace bart {

namespace system {

namespace solution {

/*! \brief Default implementation of the Solution class.
 *
 * Although this implements the methods of SolutionI it does not have any
 * specified template values so it continues to be abstract.
 *
 * @tparam IndexType type of index used to identify solutions
 * @tparam SolutionType type for solutions.
 */
template <typename IndexType, typename SolutionType>
class Solution : public SolutionI<IndexType, SolutionType> {
 public:
  using typename SolutionI<IndexType, SolutionType>::SolutionMap;
  virtual ~Solution() = default;

  virtual const SolutionMap& solutions() const override {
    return solutions_;
  }

  virtual const SolutionType& operator[](const IndexType index) const override {
    return solutions_.at(index);
  }

  virtual       SolutionType& operator[](const IndexType index) override {
    return solutions_.at(index);
  }

 protected:
  SolutionMap solutions_ = {};
};

} // namespace solution

} // namespace system

} // namespace bart

#endif //BART_SRC_SYSTEM_SOLUTION_SOLUTION_H_
