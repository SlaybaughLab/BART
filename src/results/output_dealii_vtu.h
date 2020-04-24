#ifndef BART_SRC_RESULTS_OUTPUT_DEALII_VTU_H_
#define BART_SRC_RESULTS_OUTPUT_DEALII_VTU_H_

#include <memory>

#include <deal.II/numerics/data_out.h>

#include "domain/definition_i.h"
#include "results/output_i.h"

namespace bart {

namespace results {

template <int dim>
class OutputDealiiVtu : public OutputI {
 public:
  OutputDealiiVtu(const std::shared_ptr<domain::DefinitionI<dim>> &domain_ptr);

  void AddData(system::System &to_output) override;
  void WriteData(std::ostream &output_stream) const override;
  void WriteMasterFile(std::ostream &output_stream,
                       std::vector<std::string> filenames) const override;

  domain::DefinitionI<dim>* domain_ptr() const { return domain_ptr_.get(); };

 private:
  dealii::DataOut<dim> data_out_;
  std::shared_ptr<domain::DefinitionI<dim>> domain_ptr_ = nullptr;
};

} // namespace result

} //namespace bart

#endif //BART_SRC_RESULTS_OUTPUT_DEALII_VTU_H_
