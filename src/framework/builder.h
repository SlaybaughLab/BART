#ifndef BART_SRC_FRAMEWORK_BUILDER_H_
#define BART_SRC_FRAMEWORK_BUILDER_H_

#include <memory>

#include "../problem/parameters_i.h"
#include "../domain/definition.h"

namespace bart {

namespace framework {

template <int dim>
class Builder {
 public:
  Builder(std::shared_ptr<problem::ParametersI> parameters);
  ~Builder() = default;

 private:
  std::shared_ptr<problem::ParametersI> parameters_;
  std::unique_ptr<domain::Definition<dim>> BuildDefinition() const;
  
};

} // namespace framework

} // namespace bart 

#endif // BART_SRC_FRAMEWORK_BUILDER_H_
