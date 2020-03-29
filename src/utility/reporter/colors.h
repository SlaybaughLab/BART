#ifndef BART_SRC_UTILITY_REPORTER_COLORS_H_
#define BART_SRC_UTILITY_REPORTER_COLORS_H_

#include <unordered_map>
#include <string>

namespace bart {

namespace utility {

namespace reporter {

enum class Color {
  Reset = 0,
  Red = 1,
  Green = 2,
  Blue = 3,
};

} // namespace reporter

} // namespace utility

} // namespace bart

#endif //BART_SRC_UTILITY_REPORTER_COLORS_H_
