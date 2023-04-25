#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::SparseMatrix;                  // variable size matrix, double precision

// [[Rcpp::export]]
MatrixXd sparse_mat_mult(
    const SparseMatrix<double> M_sparse,
    MatrixXd M) {
  return (M_sparse * M);
}

// [[Rcpp::export]]
MatrixXd sparse_mat_chol(
    const SparseMatrix<double> M_sparse,
    MatrixXd M) {
  return (M_sparse * M);
}