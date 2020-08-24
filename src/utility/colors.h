#ifndef BART_SRC_UTILITY_COLORS_H_
#define BART_SRC_UTILITY_COLORS_H_

#include <map>
#include <string>

namespace bart {

namespace utility {

enum class Color {
  kReset = 0,
  kRed = 1,
  kGreen = 2,
  kBlue = 3,
  kYellow = 4,
};

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

} // namespace utility

} // namespace bart

#endif //BART_SRC_UTILITY_COLORS_H_
