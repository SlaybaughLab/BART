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

#define BART_INTERFACE_FACTORY(interface, names) \
template <typename ...T> \
class interface##Factory : public utility::factory::AutoRegisteringFactory<names, \
std::unique_ptr<interface>(*)(T...)> {};

#define BART_INTERFACE_FACTORY_REGISTRAR(interface, names) \
template <typename T, typename ...U> \
class interface##FactoryRegistrar { \
public: \
using interface##Constructor = std::unique_ptr<interface>(*)(U...); \
interface##FactoryRegistrar(const names name, \
                           const interface##Constructor& constructor) { \
    interface##Factory<U...>::get().RegisterConstructor(name, constructor); \
} \
};

//#define BART_REGISTER_CONSTRUCTOR(interface, class_name, )


} // namespace factory

} // namespace utility

} // namespace bart

#endif //BART_SRC_UTILITY_FACTORY_AUTO_REGISTERING_FACTORY_H_
