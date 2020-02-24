#ifndef BART_SRC_FORMULATION_CFEM_DIFFUSION_STAMPER_H_
#define BART_SRC_FORMULATION_CFEM_DIFFUSION_STAMPER_H_

#include <functional>
#include <memory>
#include <unordered_set>
#include <vector>

#include "system/moments/spherical_harmonic_types.h"
#include "domain/definition_i.h"
#include "formulation/scalar/cfem_diffusion_i.h"
#include "formulation/cfem_stamper_i.h"
#include "problem/parameter_types.h"

namespace bart {

namespace formulation {

template <int dim>
class CFEM_DiffusionStamper : public CFEMStamperI {
 public:
  using GroupNumber = int;
  using MPISparseMatrix = dealii::PETScWrappers::MPI::SparseMatrix;
  using MPIVector = dealii::PETScWrappers::MPI::Vector;
  using Boundary = bart::problem::Boundary;
  using BoundaryType = typename formulation::scalar::CFEM_DiffusionI<dim>::BoundaryType;

  /*! Constructor, using a provided set of reflective boundary conditions.
   *
   * Maintains ownership of a problem domain, and some CFEM Diffusion
   * formulation.
   *
   * @param diffusion_ptr pointer to diffusion formulation.
   * @param definition_ptr pointer to problem definition.
   * @param reflective_boundaries unordered set containing system reflective
   *        boundaries.
   */
  CFEM_DiffusionStamper(
      std::unique_ptr<formulation::scalar::CFEM_DiffusionI<dim>> diffusion_ptr,
      std::shared_ptr<domain::DefinitionI<dim>> definition_ptr,
      const std::unordered_set<Boundary> reflective_boundaries = {});

  /*! Constructor, using a provided set of reflective boundary conditions.
   *
   * This constructor accepts a mapping of each boundary to a bool indicating
   * if it is reflective. This is compatible with the output of
   * bart::problem::ParametersI::ReflectiveBoundary().
   *
   * @param diffusion_ptr pointer to diffusion formulation.
   * @param definition_ptr pointer to problem definition.
   * @param reflective_boundary_map mapping of reflective boundaries.
   */
  CFEM_DiffusionStamper(
      std::unique_ptr<formulation::scalar::CFEM_DiffusionI<dim>> diffusion_ptr,
      std::shared_ptr<domain::DefinitionI<dim>> definition_ptr,
      const std::map<Boundary, bool> reflective_boundary_map);


  void StampStreamingTerm(MPISparseMatrix& to_stamp,
                          const GroupNumber group) override;
  void StampCollisionTerm(MPISparseMatrix& to_stamp,
                          const GroupNumber group) override;
  void StampBoundaryTerm(MPISparseMatrix& to_stamp) override;
  void StampFixedSource(MPIVector& to_stamp, const GroupNumber group) override;
  void StampFissionSource(MPIVector& to_stamp,
                          const GroupNumber group,
                          const double k_effective,
                          const system::moments::MomentVector& in_group_moment,
                          const system::moments::MomentsMap& group_moments) override;
  void StampScatteringSource(MPIVector &to_stamp,
                             const GroupNumber group,
                             const system::moments::MomentsMap &group_moments) override;

  CFEM_DiffusionStamper& AddReflectiveBoundary(Boundary boundary) override {
    reflective_boundaries_.insert(boundary);
    return *this;
  }

  CFEM_DiffusionStamper& RemoveReflectiveBoundary(Boundary boundary)  override {
    reflective_boundaries_.erase(boundary);
    return *this;
  }


  std::unordered_set<Boundary> reflective_boundaries() const override {
    return reflective_boundaries_; };

 private:
  using InitializationToken =
      typename formulation::scalar::CFEM_DiffusionI<dim>::InitializationToken;

  void StampMatrix(
      MPISparseMatrix& to_stamp,
      std::function<void(dealii::FullMatrix<double>&,
                         const domain::CellPtr<dim>&)> function);

  void StampVector(
      MPIVector& to_stamp,
      std::function<void(dealii::Vector<double>&,
                         const domain::CellPtr<dim>&)> function);

  std::unique_ptr<formulation::scalar::CFEM_DiffusionI<dim>> diffusion_ptr_;
  std::shared_ptr<domain::DefinitionI<dim>> definition_ptr_;
  InitializationToken diffusion_init_token_;
  std::vector<domain::CellPtr<dim>> cells_;
  std::unordered_set<Boundary> reflective_boundaries_ = {};
};

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_CFEM_DIFFUSION_STAMPER_H_