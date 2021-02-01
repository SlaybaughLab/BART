#ifndef BART_SRC_INSTRUMENTATION_OUTSTREAM_VECTOR_TO_VTU_HPP_
#define BART_SRC_INSTRUMENTATION_OUTSTREAM_VECTOR_TO_VTU_HPP_

#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>

#include "domain/definition_i.h"
#include "instrumentation/outstream/outstream_i.h"

namespace bart::instrumentation::outstream {

template <int dim>
class VectorToVTU : public OutstreamI<dealii::Vector<double>> {
 public:
  using Vector = dealii::Vector<double>;
  using Definition = domain::DefinitionI<dim>;
  VectorToVTU(std::shared_ptr<Definition>, std::string data_name, std::string directory, std::string filename_base);
  auto Output(const Vector& to_output) -> VectorToVTU& override;

  auto definition_ptr() -> Definition* { return definition_ptr_.get(); }

  auto data_name() -> std::string { return data_name_; }
  auto directory() -> std::string { return directory_; }
  auto filename_base() -> std::string { return filename_base_; }
 private:
  std::shared_ptr<Definition> definition_ptr_{ nullptr };
  const std::string data_name_;
  const std::string directory_;
  const std::string filename_base_;

  int counter_{ 0 };
};

} // namespace bart::instrumentation::outstream

#endif //BART_SRC_INSTRUMENTATION_OUTSTREAM_VECTOR_TO_VTU_HPP_
