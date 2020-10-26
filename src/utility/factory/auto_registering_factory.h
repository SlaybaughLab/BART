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

  static auto get() -> AutoRegisteringFactory& {
    static AutoRegisteringFactory instance;
    return instance;
  }

  static auto RegisterConstructor(const ClassName name, const Constructor& constructor) -> bool {
    return constructors().insert(std::make_pair(name, constructor)).second;
  }

  static auto GetConstructor(const ClassName name) -> Constructor {
    try {
      return constructors().at(name);
    } catch (std::out_of_range&) {
      throw(std::out_of_range("Error retrieving constructor"));
    }
  }

  static auto constructors() -> std::unordered_map<ClassName, Constructor>& {
    static std::unordered_map<ClassName, Constructor> constructors_;
    return constructors_;
  }

 private:
  AutoRegisteringFactory() = default;
  AutoRegisteringFactory(const AutoRegisteringFactory&) = default;
};

#define BART_INTERFACE_FACTORY(interface, names) \
template <typename ...T> \
class interface##Factory : public utility::factory::AutoRegisteringFactory<names, \
std::unique_ptr<interface>(*)(T...)> {};

} // namespace factory

} // namespace utility

} // namespace bart

#endif //BART_SRC_UTILITY_FACTORY_AUTO_REGISTERING_FACTORY_H_
