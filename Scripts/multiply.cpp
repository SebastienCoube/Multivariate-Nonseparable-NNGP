
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]

//  // [[Rcpp::export]]
// NumericVector Multiply_vector_compressed_symmat(NumericVector v, NumericVector compressed_symmat, IntegerMatrix idx_mat, int n) 
// {
//   NumericVector res(n);
//   for(int i=0; i<compressed_symmat.length(); i ++)
//   {
//     res(idx_mat(i, 1)-1) = res(idx_mat(i, 1)-1) + v(idx_mat(i, 0)-1) * compressed_symmat(i);
//     if(idx_mat(i, 0)-1 != idx_mat(i, 1)-1)
//     {
//       res(idx_mat(i, 0)-1) = res(idx_mat(i, 0)-1) + v(idx_mat(i, 1)-1) * compressed_symmat(i);
//     }
//   }
//   return res;
// }

// [[Rcpp::export]]
NumericVector Multiply_vector_compressed_symmat(NumericVector v, NumericVector compressed_symmat, IntegerMatrix idx_mat, int n) 
{
  NumericVector res(n);
  for(int i=0; i<idx_mat.nrow(); i ++)
  {
    int j = idx_mat(i, 0);
    int k = idx_mat(i, 1);
    int l = (k-1)*n + j - k*(k-1)/2;
    res(k-1) = res(k-1) + v(j-1) * compressed_symmat(l-1);
    if(j-1 != k-1)
    {
      res(j-1) = res(j-1) + v(k-1) * compressed_symmat(l-1);
    }
  }
  return res;
}


// [[Rcpp::export]]
NumericMatrix Multiply_matrix_compressed_symmat(NumericMatrix v_, NumericVector compressed_symmat, IntegerMatrix idx_mat, int n) 
{
  NumericMatrix res(v_.nrow(), v_.ncol());
  for(int i=0; i<idx_mat.nrow(); i ++)
  {
    int j = idx_mat(i, 0);
    int k = idx_mat(i, 1);
    int l = (k-1)*n + j - k*(k-1)/2;
    res.row(k-1) = res.row(k-1) + v_.row(j-1) * compressed_symmat(l-1);
    if(j-1 != k-1)
    {
      res.row(j-1) = res.row(j-1) + v_.row(k-1) * compressed_symmat(l-1);
    }
  }
  return res;
}


// [[Rcpp::export]]
NumericVector Multiply_vector_side_blocks_upper_right(NumericVector v, arma::cube side_block_rectangles, IntegerMatrix idx_mat) 
{
  NumericVector res(side_block_rectangles.n_rows);
  for(int i=0; i<idx_mat.nrow(); i ++)
  {
    int j = idx_mat(i, 0);
    int k = idx_mat(i, 1);
    for(int l = 0; l< side_block_rectangles.n_slices; l++)
    {
      res(j-1) += side_block_rectangles(j-1, k-1, l) * 
        v(k-1+ l*side_block_rectangles.n_cols)
      ;
    }
  }
  return(res);
}



// [[Rcpp::export]]
NumericVector Multiply_vector_side_blocks_vertical_left(NumericVector v, arma::cube side_block_rectangles, IntegerMatrix idx_mat) 
{
  NumericVector res(side_block_rectangles.n_cols * side_block_rectangles.n_slices);
  for(int i=0; i<idx_mat.nrow(); i ++)
  {
    int j = idx_mat(i, 0);
    int k = idx_mat(i, 1);
    for(int l = 0; l< side_block_rectangles.n_slices; l++)
    {
      res(k-1 + l*side_block_rectangles.n_cols) += side_block_rectangles(j-1, k-1, l) * 
        v(j-1)
      ;
    }
  }
  return(res);
}


