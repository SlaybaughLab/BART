#ifndef BART_SRC_FORMULATION_CFEM_DIFFUSION_STAMPER_H_
#define BART_SRC_FORMULATION_CFEM_DIFFUSION_STAMPER_H_

#include <functional>
#include <memory>
#include <vector>

#include "domain/definition_i.h"
#include "formulation/scalar/cfem_diffusion_i.h"
#include "formulation/stamper_i.h"
#include "problem/parameter_types.h"

namespace bart {

namespace formulation {

template <int dim>
class CFEM_DiffusionStamper : public StamperI<dim> {
 public:
  using GroupNumber = int;
  using MPISparseMatrix = dealii::PETScWrappers::MPI::SparseMatrix;
  using Boundary = bart::problem::Boundary;
  using BoundaryType = typename formulation::scalar::CFEM_DiffusionI<dim>::BoundaryType;

  CFEM_DiffusionStamper(
      std::unique_ptr<formulation::scalar::CFEM_DiffusionI<dim>> diffusion_ptr,
      std::unique_ptr<domain::DefinitionI<dim>> definition_ptr);

  void StampStreamingTerm(MPISparseMatrix& to_stamp, const GroupNumber group);
  void StampCollisionTerm(MPISparseMatrix& to_stamp, const GroupNumber group);
  void StampBoundaryTerm(MPISparseMatrix& to_stamp);

  std::vector<Boundary> reflective_boundaries() const {
    return reflective_boundaries_; };

 private:
  using InitializationToken =
      typename formulation::scalar::CFEM_DiffusionI<dim>::InitializationToken;
  using Cell = typename domain::DefinitionI<dim>::Cell;

  void StampMatrix(
      MPISparseMatrix& to_stamp,
      std::function<void(dealii::FullMatrix<double>&, const Cell&)> function);

  std::unique_ptr<formulation::scalar::CFEM_DiffusionI<dim>> diffusion_ptr_;
  std::unique_ptr<domain::DefinitionI<dim>> definition_ptr_;
  InitializationToken diffusion_init_token_;
  std::vector<Cell> cells_;
  std::vector<Boundary> reflective_boundaries_ = {};
};

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_CFEM_DIFFUSION_STAMPER_H_