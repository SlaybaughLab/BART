#ifndef BART_SRC_UTILITY_COLORS_HPP_
#define BART_SRC_UTILITY_COLORS_HPP_

#include <map>
#include <string>

namespace bart::utility {

enum class Color {
  kReset = 0,
  kRed = 1,
  kGreen = 2,
  kBlue = 3,
  kYellow = 4,
};

// LCOV_EXCL_START

inline std::string to_string(Color to_convert) {
  std::map<Color, std::string> color_string{
      {Color::kReset, "\033[0m"},
      {Color::kRed,   "\033[31m"},
      {Color::kGreen, "\033[32m"},
      {Color::kYellow,"\033[33m"},
      {Color::kBlue,  "\033[34m"},
  };
  return color_string.at(to_convert);
}

// LCOV_EXCL_STOP

} // namespace bart::utility

#endif //BART_SRC_UTILITY_COLORS_HPP_
