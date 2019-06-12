#ifndef BART_SRC_SYSTEM_SOLUTION_I_H_
#define BART_SRC_SYSTEM_SOLUTION_I_H_

#include <map>

namespace bart {

namespace system {

namespace solution {

/*! \brief Interface for solution storage class.
 *
 * This class is the interface for classes that store system solution.
 * Generally, this is just a resource management class, to own the solutions.
 * They can be defined using any time to uniquely identify them that can be
 * hashed by a mapping. Solution type is generally a vector.
 *
 * @tparam IndexType type of index used to identify solutions.
 * @tparam SolutionType type for solutions.
 */
template <typename IndexType, typename SolutionType>
class SolutionI {
 public:
  using SolutionMap = std::map<IndexType, SolutionType>;
  virtual ~SolutionI() = default;
  /*! \brief Returns the full mapping of solutions.
   *
   * This method is intended to expose the internally managed resource to allow
   * the user to use std::map methods directly.
   *
   * @return reference to the internal solutions map.
   */
  virtual const SolutionMap& solutions() const = 0;

  /*! \brief Returns the solution that corresponds to a provided index.
   *
   * @return a constant reference to the solution.
   */
  virtual const SolutionType& operator[](const IndexType) const = 0;

  /*! \brief Returns the solution that corresponds to a provided index.
   *
   * @return a reference to the solution.
   */
  virtual       SolutionType& operator[](const IndexType) = 0;
};

} // namespace solution

} // namespace system

} // namespace bart

#endif // BART_SRC_SYSTEM_SOLUTION_I_H_