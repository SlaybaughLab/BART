#ifndef BART_SRC_UTILITY_FACTORY_AUTO_REGISTERING_FACTORY_H_
#define BART_SRC_UTILITY_FACTORY_AUTO_REGISTERING_FACTORY_H_

#include <memory>
#include <unordered_map>

namespace bart {

namespace utility {

namespace factory {

template <typename ClassName, typename Constructor>
class AutoRegisteringFactory {
 public:
  //! Virtual d'tor, we expect this class to be inherited from.
  ~AutoRegisteringFactory() = default;

  static AutoRegisteringFactory& get() {
    static AutoRegisteringFactory instance;
    return instance;
  }

  bool RegisterConstructor(const ClassName name, const Constructor& constructor) {
    return constructors_.insert(std::make_pair(name, constructor)).second;
  }

  Constructor GetConstructor(const ClassName name) {
    return constructors_.at(name);
  }

 private:
  AutoRegisteringFactory() = default;
  AutoRegisteringFactory(const AutoRegisteringFactory&) = default;
  std::unordered_map<ClassName, Constructor> constructors_{};
};

} // namespace factory

} // namespace utility

} // namespace bart

#endif //BART_SRC_UTILITY_FACTORY_AUTO_REGISTERING_FACTORY_H_
