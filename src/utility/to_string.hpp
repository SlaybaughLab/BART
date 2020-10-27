#ifndef BART_SRC_UTILITY_TO_STRING_HPP_
#define BART_SRC_UTILITY_TO_STRING_HPP_

#include <string>

namespace bart::utility {

template <typename T>
inline auto to_string(T to_convert) {
  if constexpr (requires { std::to_string(to_convert); }) {
    return std::to_string(to_convert);
  } else if constexpr (requires { to_string(to_convert); }) {
    return to_string(to_convert);
  }
}

} // namespace bart::utility

#endif //BART_SRC_UTILITY_TO_STRING_HPP_
