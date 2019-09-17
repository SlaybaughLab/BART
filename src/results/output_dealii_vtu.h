#ifndef BART_SRC_RESULTS_OUTPUT_DEALII_VTU_H_
#define BART_SRC_RESULTS_OUTPUT_DEALII_VTU_H_

#include "output_i.h"

#include <memory>

#include "domain/definition_i.h"

namespace bart {

namespace results {

template <int dim>
class OutputDealiiVtu : public OutputI {
 public:
  OutputDealiiVtu(const std::shared_ptr<domain::DefinitionI<dim>> &domain_ptr)
  : domain_ptr_(domain_ptr) {}
  void Output(system::System &to_output) const override {};

  domain::DefinitionI<dim>* domain_ptr() const { return domain_ptr_.get(); };

 private:
  std::shared_ptr<domain::DefinitionI<dim>> domain_ptr_ = nullptr;
};

} // namespace result

} //namespace bart

#endif //BART_SRC_RESULTS_OUTPUT_DEALII_VTU_H_
