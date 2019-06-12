#ifndef BART_SRC_FORMULATION_CFEM_STAMPER_I_H_
#define BART_SRC_FORMULATION_CFEM_STAMPER_I_H_

#include <unordered_set>

#include "system/moments/moment_types.h"
#include "formulation/stamper_i.h"
#include "problem/parameter_types.h"

namespace bart {

namespace formulation {

/*! Stamps terms of a basic CFEM formulation onto MPI matrices and vectors.
 *
 * Functions in this class take input matrices or vectors, and stamps the
 * contribution from each cell in the mesh. In general, this class will iterate
 * over all the cells, call the formulation to get a cell matrix, and then
 * stamp that onto the MPI matrix.
 *
 * The basic stamper stamps the following terms:
 *
 * Bilinear (take a matrix as an input)
 *  - Streaming
 *  - Collision
 *  - Boundary
 *
 * Linear (take a vector as an input)
 *  - Fixed source
 *  - Fission source
 *  - Scattering source
 *
 * Source terms are calculated using moments.
 */

class CFEMStamperI : public StamperI {
 public:
  using GroupNumber = int;
  using Boundary = bart::problem::Boundary;
  using MPISparseMatrix = dealii::PETScWrappers::MPI::SparseMatrix;
  using MPIVector = dealii::PETScWrappers::MPI::Vector;

  virtual ~CFEMStamperI() = default;

  /*! Stamps all cell streaming term matrices onto the system matrix.
   *
   * @param to_stamp MPI sparse matrix to be stamped onto.
   * @param group group to stamp.
   */
  virtual void StampStreamingTerm(MPISparseMatrix& to_stamp,
                                  const GroupNumber group) = 0;

  /*! Stamps all cell collision term matrices onto the system matrix.
   *
   * @param to_stamp MPI sparse matrix to be stamped onto.
   * @param group group to stamp.
   */
  virtual void StampCollisionTerm(MPISparseMatrix& to_stamp,
                                  const GroupNumber group) = 0;

  /*! Stamps all cell boundary term matrices onto the system matrix.
   *
   * @param to_stamp MPI sparse matrix to be stamped onto.
   */
  virtual void StampBoundaryTerm(MPISparseMatrix& to_stamp) = 0;

  /*! Stamps all cell fixed source terms onto a system vector.
   *
   * @param to_stamp MPI vector to be stamped onto.
   * @param group group to stamp.
   */
  virtual void StampFixedSource(MPIVector& to_stamp,
                                const GroupNumber group) = 0;

  /*! Stamps all cell fission source terms onto a system vector.
   *
   * @param to_stamp MPI vector to be stamped onto
   * @param group group to stamp.
   * @param k_effective k_effective for the system.
   * @param in_group_moment moment to be used for the in-group contribution.
   * @param group_moments moments to be used for all out-group contributions.
   */
  virtual void StampFissionSource(MPIVector& to_stamp,
                                  const GroupNumber group,
                                  const double k_effective,
                                  const data::MomentVector& in_group_moment,
                                  const data::MomentsMap& group_moments) = 0;

  /*! Stamps all cell scattering source terms onto a system vector.
   *
   * @param to_stamp MPI vector to be stamped onto
   * @param group group to stamp.
   * @param in_group_moment moment to be used for the in-group contribution.
   * @param group_moments moments to be used for all out-group contributions.
   */
  virtual void StampScatteringSource(MPIVector& to_stamp,
                                     const GroupNumber group,
                                     const data::MomentVector& in_group_moment,
                                     const data::MomentsMap& group_moments) = 0;


  /*! Adds a system reflective boundary.
   *
   * Note that you can add boundaries for dimensions that aren't appropriate
   * for the problem (e.g. in the z-direction for a 1D problem) with no effect.
   *
   * @param boundary boundary to add.
   * @return this object
   */
  virtual CFEMStamperI& AddReflectiveBoundary(Boundary boundary) = 0;

  /*! Removes a system reflective boundary.
   *
   * @param boundary boundary to remove.
   * @return this object
   */
  virtual CFEMStamperI& RemoveReflectiveBoundary(Boundary boundary) = 0;

  /*! Returns the reflective boundary conditions.
   *
   * @return an unordered set holding the reflective boundaries.
   */
  virtual std::unordered_set<Boundary> reflective_boundaries() const = 0;

};

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_CFEM_STAMPER_I_H_