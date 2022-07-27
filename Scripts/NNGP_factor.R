
NNGP_factor  = function(locs, NNarray, covparms, covfun){
  rows = apply(NNarray, 1, function(NNarray_row){
    NNarray_row = c(na.omit(NNarray_row))
    return(solve(chol(covfun(covparms, locs[NNarray_row[seq(length(NNarray_row), 1)],, drop = F])))[seq(length(NNarray_row), 1), length(NNarray_row)])
  })
  Linv = matrix(NA, length(rows), max(sapply(rows, length)))
  Linv[cbind(rep(seq(length(rows)), sapply(rows, length)), unlist(lapply(rows, function(x)seq(length(x)))))] = unlist(rows)
  return(Linv)
}



locs = matrix(runif(10000), ncol = 2)
NNarray = GpGp::find_ordered_nn(locs, 10)
Linv = NNGP_factor(locs, NNarray, c(1, .1, 0), GpGp::matern15_isotropic)
L  = Matrix::sparseMatrix(i = row(NNarray)[!is.na(NNarray)], j = NNarray[!is.na(NNarray)], x = Linv[!is.na(NNarray)])
Bidart::plot_pointillist_painting(locs, as.vector(Matrix::solve(L, rnorm(nrow(locs)))))



covfun = GpGp::matern15_isotropic
covparms = c(1, .1, 0)

NNGP_factor_and_diff  = function(locs, NNarray, covparms, covfun){
  for(i in seq(nrow(NNarray))){
    locsub = na.omit(NNarray[i,])
    locsub = locsub[seq(length(locsub), 1)]
    cov_chol = chol(covfun(covparms, locs[locsub,,drop = F]))
  }
}