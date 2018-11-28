#include "material_base.h"

static MaterialBase::dealii::Vector<double> GetEigenVectors (
    const dealii::FullMatrix<double>& mat) const {
  AssertThrow (mat.m()==mat.n(), ExcMessage("matrix must be squared"));
  // the rest is using power iteration to calculate eigenvectors. The assumption
  // is that we definitely have a dominant mode.

  // we are solving: Ab = lambda*b
  // lambda = b^T * A * b / b^T * b
  // b = A*b / ||A*b||
  dealii::Vector<double> v1(mat.n()),v2;
  // const dealii::FullMatrix<double> id_mat(IdentityMatrix(mat.n()));
  v1 = 1.;
  double err = 1.;
  double eig = mat.matrix_norm_square(v1) / v1.norm_sqr();
  while (err>1.e-10) {
    double eig1 = eig;
    // b = A*b / ||A*b||
    mat.vmult(v2, v1);// note that v2 is required to be different from v1 by dealii
    auto inv_scal = 1./v2.l1_norm();
    v2 *= inv_scal;
    v1 = v2;
    // lambda = b^T * A * b / b^T * b
    eig = mat.matrix_norm_square(v1) / v1.norm_sqr();
    // calculate error between adjacent iterations
    err = std::fabs(eig1 - eig) / std::fabs(eig);
  }
  return v1;
}
