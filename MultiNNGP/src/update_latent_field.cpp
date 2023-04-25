#include <RcppEigen.h>
#include <Eigen/Sparse>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::MatrixXd;
using Eigen::SparseMatrix;
typedef Eigen::SparseMatrix<double> SpMat;



SpMat blockToeplitz(std::vector<SpMat> blocks) {
  int num_blocks = blocks.size();
  int block_size = blocks[0].rows();
  int matrix_size = num_blocks * block_size;
  
  std::vector<Eigen::Triplet<double>> tripletList;
  tripletList.reserve(matrix_size * matrix_size);
  
  for (int i = 0; i < matrix_size; i += block_size) {
    for (int j = 0; j < matrix_size; j += block_size) {
      int b = i / block_size; // block index
      int k = i - b * block_size; // row index within block
      
      for (int p = 0; p < block_size; ++p) {
        int q = (j + p - k) % block_size; // column index within block
        int j_prime = j + q;
        int b_prime = j_prime / block_size; // block index
        tripletList.emplace_back(i + p, j_prime, blocks[(b + b_prime) % num_blocks].coeff(p, q));
      }
    }
  }
  
  SpMat result(matrix_size, matrix_size);
  result.setFromTriplets(tripletList.begin(), tripletList.end());
  
  return result;
}



// [[Rcpp::export]]
int update_latent_field_markov_inputs(
    Rcpp::List noise_precisions, 
    MatrixXd y_minus_field
) {
  return (1);
}