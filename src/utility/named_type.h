#ifndef BART_SRC_UTILITY_NAMED_TYPE_H_
#define BART_SRC_UTILITY_NAMED_TYPE_H_

#include <memory>

namespace bart {

namespace utility {

/*! \brief Class providing strong types for generic types.
 *
 * Drawn directly from Jonathan Boccara's blog, <a href="https://www.fluentcpp.com/2016/12/08/strong-types-for-strong-interfaces/">Fluent C++</a>.
 *
 * \tparam T type to be wrapped.
 * \tparam Parameter a parameter unique to the type, to prevent types just
 * being aliases of each other.
 *
 * See the website for more specific usage.
 */
template <typename T, typename Parameter>
class NamedType
{
 public:
  explicit NamedType(T const& value) : value_(value) {}
  explicit NamedType(T&& value) : value_(std::move(value)) {}
  T& get() { return value_; }
  T const& get() const {return value_; }
  bool operator==(const NamedType& rhs) const {
    return rhs.get() == value_;
  }
  bool operator<(const NamedType& rhs) const {
    return rhs.get() < value_;
  }
  template <typename U = T,
      std::enable_if_t<std::is_same_v<U, bool>, nullptr_t> = nullptr>
  operator bool() const { return value_; }
 private:
  T value_;
};

} // namespace utility

} // namespace bart

#endif // BART_SRC_UTILITY_NAMED_TYPE_H_