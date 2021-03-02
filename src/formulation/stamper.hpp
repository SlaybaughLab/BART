#ifndef BART_SRC_FORMULATION_STAMPER_HPP_
#define BART_SRC_FORMULATION_STAMPER_HPP_

#include <memory>

#include "domain/domain_i.hpp"
#include "formulation/stamper_i.hpp"
#include "utility/has_dependencies.h"

namespace bart::formulation {

/*! \brief Default implementation of the domain stamper.
 *
 *  \author J.S. Rehak
 *  \tparam dim spatial dimension of the cells in the mesh
 */
template <int dim>
 class Stamper : public StamperI<dim>, public utility::HasDependencies {
 public:
  using typename StamperI<dim>::CellMatrixStampFunction;
  using typename StamperI<dim>::CellVectorStampFunction;
  using typename StamperI<dim>::FaceMatrixStampFunction;
  using typename StamperI<dim>::FaceVectorStampFunction;
  using Domain = typename domain::DomainI<dim>;
  /*! \brief Constructor.
   * Takes a domain definition dependency that provides the domain of cells to iterate over. The matrices and vectors
   * passed in this classes functions should be seperately initialized using this domain. */
  explicit Stamper(std::shared_ptr<Domain>);
  virtual ~Stamper() = default;

  auto StampMatrix(system::MPISparseMatrix& to_stamp, CellMatrixStampFunction stamp_function) -> void override;
  auto StampVector(system::MPIVector& to_stamp, CellVectorStampFunction stamp_function) -> void override;
  auto StampBoundaryMatrix(system::MPISparseMatrix &to_stamp, FaceMatrixStampFunction stamp_function) -> void override;
  auto StampBoundaryVector(system::MPIVector &to_stamp, FaceVectorStampFunction stamp_function) -> void override;

  /*! \brief Access domain definition dependency */
  auto domain_ptr() const { return domain_ptr_.get(); }
 private:
  std::shared_ptr<Domain> domain_ptr_;
};

} // namespace bart::formulation

#endif //BART_SRC_FORMULATION_STAMPER_HPP_
