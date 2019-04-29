#ifndef BART_SRC_FORMULATION_CFEM_DIFFUSION_STAMPER_H_
#define BART_SRC_FORMULATION_CFEM_DIFFUSION_STAMPER_H_

#include <functional>
#include <memory>
#include <unordered_set>
#include <vector>

#include "data/moment_types.h"
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
  using MPIVector = dealii::PETScWrappers::MPI::Vector;
  using Boundary = bart::problem::Boundary;
  using BoundaryType = typename formulation::scalar::CFEM_DiffusionI<dim>::BoundaryType;

  CFEM_DiffusionStamper(
      std::unique_ptr<formulation::scalar::CFEM_DiffusionI<dim>> diffusion_ptr,
      std::unique_ptr<domain::DefinitionI<dim>> definition_ptr,
      const std::unordered_set<Boundary> reflective_boundaries = {});
  CFEM_DiffusionStamper(
      std::unique_ptr<formulation::scalar::CFEM_DiffusionI<dim>> diffusion_ptr,
      std::unique_ptr<domain::DefinitionI<dim>> definition_ptr,
      const std::map<Boundary, bool> reflective_boundary_map);


  void StampStreamingTerm(MPISparseMatrix& to_stamp, const GroupNumber group);
  void StampCollisionTerm(MPISparseMatrix& to_stamp, const GroupNumber group);
  void StampBoundaryTerm(MPISparseMatrix& to_stamp);
  void StampFixedSource(MPIVector& to_stamp, const GroupNumber group);
  void StampFissionSource(MPIVector& to_stamp,
                          const GroupNumber group,
                          const double k_effective,
                          const data::MomentVector& in_group_moment,
                          const data::MomentsMap& group_moments);

  CFEM_DiffusionStamper& AddReflectiveBoundary(Boundary boundary) {
    reflective_boundaries_.insert(boundary);
    return *this;
  }

  CFEM_DiffusionStamper& RemoveReflectiveBoundary(Boundary boundary) {
    reflective_boundaries_.erase(boundary);
    return *this;
  }


  std::unordered_set<Boundary> reflective_boundaries() const {
    return reflective_boundaries_; };

 private:
  using InitializationToken =
      typename formulation::scalar::CFEM_DiffusionI<dim>::InitializationToken;
  using Cell = typename domain::DefinitionI<dim>::Cell;

  void StampMatrix(
      MPISparseMatrix& to_stamp,
      std::function<void(dealii::FullMatrix<double>&, const Cell&)> function);

  void StampVector(
      MPIVector& to_stamp,
      std::function<void(dealii::Vector<double>&, const Cell&)> function);

  std::unique_ptr<formulation::scalar::CFEM_DiffusionI<dim>> diffusion_ptr_;
  std::unique_ptr<domain::DefinitionI<dim>> definition_ptr_;
  InitializationToken diffusion_init_token_;
  std::vector<Cell> cells_;
  std::unordered_set<Boundary> reflective_boundaries_ = {};
};

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_CFEM_DIFFUSION_STAMPER_H_