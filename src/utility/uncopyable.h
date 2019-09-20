#ifndef BART_SRC_UTILITY_UNCOPYABLE_H_
#define BART_SRC_UTILITY_UNCOPYABLE_H_

/*! A base class that declares the copy and copy-assignment constructors private
 * to prevent unneeded or bad copies. */

namespace bart {

namespace utility {

class Uncopyable {
 public:
  Uncopyable(const Uncopyable&) = delete;
  Uncopyable& operator=(const Uncopyable&) = delete;
 protected:
  Uncopyable() = default;
  ~Uncopyable() = default;
};

} // namespace utility

} // namespace bart
  
 

#endif // BART_SRC_UTILITY_UNCOPYABLE_H_
