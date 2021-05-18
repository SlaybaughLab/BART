#ifndef BART_SRC_INSTRUMENTATION_OUTSTREAM_VECTOR_MAP_TO_VTU_HPP_
#define BART_SRC_INSTRUMENTATION_OUTSTREAM_VECTOR_MAP_TO_VTU_HPP_

#include <memory>
#include <unordered_map>

#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>

#include "domain/domain_i.hpp"
#include "instrumentation/outstream/outstream_i.h"
#include "instrumentation/outstream/factory.hpp"

namespace bart::instrumentation::outstream {

template <int dim>
class VectorMapToVTU : public OutstreamI<std::unordered_map<int, dealii::Vector<double>>> {
 public:
  using Vector = dealii::Vector<double>;
  using Definition = domain::DomainI<dim>;
  using VectorMap = std::unordered_map<int, Vector>;

  VectorMapToVTU(std::shared_ptr<Definition>, std::string data_name, std::string directory, std::string filename_base);
  auto Output(const VectorMap& to_output) -> VectorMapToVTU& override;

  auto definition_ptr() const { return definition_ptr_.get(); }
  auto data_name() const { return data_name_; }
  auto directory() const { return directory_; }
  auto filename_base() const { return filename_base_; }
 private:
  std::shared_ptr<Definition> definition_ptr_{ nullptr };
  const std::string data_name_;
  const std::string directory_;
  const std::string filename_base_;

  int counter_{ 0 };
  static bool is_registered_;
};

} // namespace bart::instrumentation::outstream

#endif //BART_SRC_INSTRUMENTATION_OUTSTREAM_VECTOR_MAP_TO_VTU_HPP_
