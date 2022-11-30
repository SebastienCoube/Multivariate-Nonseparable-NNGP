#################
# Generate data #
#################

# generate test data
set.seed(1)
n_var = 10


locs_1 = 10*cbind(runif(40), runif(40))
locs_2 = 10*cbind(runif(50), runif(50))
var_tag_1 = 1+floor(n_var*runif(40))
var_tag_2 = 1+floor(n_var*runif(50))

rho = GpGp::exponential_isotropic(c(1, 1, 0), .1*matrix(rnorm(2*n_var), n_var))
rho_vec = rho[lower.tri(rho, diag = F)];  remove(rho)

a2_vec = runif(n_var)
nu_vec = .5 + 2*runif(n_var)
alpha = .0001 * runif(1)
a = runif(1) 
b = runif(1) 
cc = .1 * runif(1) 
lambda = runif(1) 
delta = runif(1) 
r = .1 * runif(1) 
A_vec = runif(n_var)
u = seq(0, 5)


######################################
# functions to get components of GMA #
######################################

get_multiplier = function(
    a, b, cc, delta, lambda, r, 
    A_vec, nu_vec, a2_vec, 
    u
){
  A_ = outer(A_vec, A_vec, "*") ; A_ = A_[lower.tri(A_, T)]
  a_ = outer(a2_vec, a2_vec, "+") ; a_ =a_/2; a_ = a_[lower.tri(a_, T)]; a_ = sqrt(a_)
  nu_ = outer(nu_vec, nu_vec, "+") ; nu_ = nu_/2; nu_ = nu_[lower.tri(nu_, T)]
  tau_ = a2_vec^(nu_vec)*(1-A_vec^2)^2 / gamma(nu_vec) # n_var 18 τii = σi2 a2νii (1 − Ai ) (here sigma is part of rho, added later)
  tau_ = outer(sqrt(tau_), sqrt(tau_)); tau_ = tau_[lower.tri(tau_, T)]
  res=  
    (
      outer(
        rep(1, length(u)), 
        (2^(1-nu_)) * 
          tau_
        * #gamma(nu_)*
          a_^(-2*nu_)
      )
    ) / (
      (
        outer((1 + (cc*u)^(2*a)     )^  delta,  rep(1, length(A_))) - 
          outer((1 + (r *u)^(2*lambda))^(-delta), A_                )
      )
      * 
        (
          outer((1 + (cc*u)^(2*a)     )^  b,  rep(1, length(A_))) - 
            outer((1 + (r *u)^(2*lambda))^(-b), A_                )
        )
    )
  row.names(res) = paste("u=", u, sep="")
  M = matrix(0, length(A_vec),length(A_vec))
  colnames(res) = paste("var", row(M)[lower.tri(M, T)], col(M)[lower.tri(M, T)] , sep="_")
  res
}

get_effective_range = function(
    a, b, cc, delta, lambda, r, 
    A_vec, a2_vec, 
    u
)
{
  A_ = outer(A_vec, A_vec, "*") ; A_ = A_[lower.tri(A_, T)]
  a_ = outer(a2_vec, a2_vec, "+") ; a_ =a_/2; a_ = a_[lower.tri(a_, T)]; a_ = sqrt(a_)
  outer(
    rep(1, length(u)), 
    a_
  )/
    (
      outer((1 + (cc*u)^(2*a)     )^  b,  rep(1, length(A_))) - 
        outer((1 + (r *u)^(2*lambda))^(-b), A_                )
    )^.5
}

Matern = function(h, r, nu_)(r * h)^nu_ * besselK(r * h, nu_)

##################################################
# Some functions to work with symmetric matrices #
##################################################

# given integer n_loc, gives 2-column matrix of lower triangular indices 
# in the square matrix n_loc*n_loc WITHOUT diagonal as in dist(,diag = F)


# ex :              index | i j
# n_loc= 4 -> 1 . . .    1    | 1 1
#         2 5 . .    2    | 2 1
#         3 6 8 .    3    | 3 1
#         4 7 9 10   4    | 4 1
#                    5    | 2 2
#                    6    | 3 2 
#                    7    | 4 2 
#                    8    | 3 3 
#                    9    | 4 3
#                    10   | 4 4 
#                      

lower_tri_idx = function(n_loc, diag = F)
{
  
  if(diag == F)
  {
    return(
      cbind(
        rev(abs(sequence(seq.int(n_loc - 1)) - n_loc) + 1),
        rep.int(seq.int(n_loc - 1), rev(seq.int(n_loc - 1)))
      )
    )
  }
  if(diag == T)
  {
    return(
      cbind(
        rev(abs(sequence(seq.int(n_loc)) - n_loc) + 1),
        rep.int(seq.int(n_loc), rev(seq.int(n_loc)))
      )
    )
  }
}
lower_tri_idx(10)

# given three integers i, j, n_loc, gives the position of coefficient (i, j)
# in the n_loc*n_loc square matrix in the lower triangular coefficients
# WITH diagonal as in dist(,diag = T)

# ex: i = 4, j = 3, n_loc = 4
#             j
#         1 . . . 
#         2 5 . . 
#         3 6 8 . 
#       i 4 7(9)10
#     result = 9

position_in_lower_tri = function(i, j, n_loc)
{  
  a = pmin(i, j); b = pmax(i, j)
  (a-1)*(n_loc-a+1) + (a-1)*a/2 + b-a+1
}
#  i = 1;j =3;n_loc=3;position_in_lower_tri(i, j, n_loc)
#  i = 2;j =1;n_loc=3;position_in_lower_tri(i, j, n_loc)
#  i = 2;j =3;n_loc=3;position_in_lower_tri(i, j, n_loc)
#  i = 3;j =3;n_loc=3;position_in_lower_tri(i, j, n_loc)


# given a vector of integers i_vec and an integer n_var, 
# n_loc var being greater than the entries of i_vec, 
# computes the lower-triangular without diagonal of the crossing matrix of i_vec, i_vec
# and finds the position of those pairs in the lower triangular part of the 
# matrix of size n_var * n_var
position_in_lower_tri_cross_vec = function(i_vec, n_var, diag = F)
{  
  idx = lower_tri_idx(length(i_vec), diag = diag)
  position_in_lower_tri(i_vec[idx[,1]], i_vec[idx[,2]], n_var)
}

##  i_vec = 1 + floor(n_var*runif(10))
##  position_in_lower_tri_cross_vec(i_vec, n_var)
##  length(position_in_lower_tri_cross_vec(i_vec, n_var))
##  
##  i_vec = 1 + floor(n_var*runif(10))
##  position_in_lower_tri_cross_vec(i_vec, n_var)
##  length(position_in_lower_tri_cross_vec(i_vec,n_var))


put_ones_in_rho_vec = function(rho_vec)
{
  nn = .5 + sqrt(.25 + 2*length(rho_vec))
  M = matrix(1, nn, nn)
  M[lower.tri(M)]=rho_vec
  M[lower.tri(M, T)]
}

######################
# getting GMA covmat #
######################

# Toy example  
multiplier = get_multiplier(
  a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
  r = r, A_vec = A_vec, nu_vec = nu_vec, a2_vec = a2_vec, u = u
)
effective_range = get_effective_range(
  a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
  r = r, A_vec = A_vec, a2_vec = a2_vec, u = u
)


# compute unique elements of GMA covmat
GMA_compressed = function(
    locs, 
    var_tag,
    multiplier, 
    effective_range,
    nu_, 
    rho_vec_with_ones, 
    n_var
)
{
  # size of spatial sets of interest
  n_loc = nrow(locs)
  n_lags = nrow(multiplier)
  # var combinations
  var_idx = position_in_lower_tri_cross_vec(var_tag, n_var, diag = T)
  # distance between locs pairs
  h = as.matrix(dist(locs, diag = T));h = h[lower.tri(h, diag = T)]
  h[h==0] = min(h[h!=0])*.0001
  # indices of pairs in covariance matrix
  idx = lower_tri_idx(n_loc, diag = T);idx[,1] = idx[,1]
  # creating matrices
  res = list()
  res$covmat_without_rho = matrix(0, nrow(idx), n_lags)
  res$covmat = matrix(0, nrow(idx), n_lags)
  
  # looping over time lags
  for(i_lag in seq(n_lags)){
    # index of position in the covariance matrix
    # covariance without correlation
    res$covmat_without_rho[,i_lag] = 
      Matern(h, r = effective_range[i_lag,var_idx], nu_ = nu_[var_idx]) * 
      multiplier[i_lag, var_idx]
    # covariance
    res$covmat[,i_lag] = res$covmat_without_rho[,i_lag] * rho_vec_with_ones[var_idx]
  }
  res
}




# elements of rectangular GMA matrix
GMA_rectangular = function(
    locs_1, locs_2, 
    var_tag_1, var_tag_2,
    multiplier, 
    effective_range,
    nu_, 
    rho_vec_with_ones, 
    n_var
)
{
  # size of spatial sets of interest
  n_loc_1 = nrow(locs_1)
  n_loc_2 = nrow(locs_1)
  n_lags = nrow(multiplier)
  # var combinations
  tatato = expand.grid(var_tag_1, var_tag_2)
  var_idx = position_in_lower_tri(tatato[,1], tatato[,2], n_var)
  # distance between locs pairs
  h = fields::rdist(locs_1, locs_2)
  h[h==0] = min(h[h!=0])*.0001
  # creating matrices
  res = list()
  res$covmat_without_rho = array(0, c(dim(h), nrow(multiplier)))
  res$covmat = array(0, c(dim(h), nrow(multiplier)))
  
  # looping over time lags
  for(i_lag in seq(n_lags)){
    # index of position in the covariance matrix
    # covariance without correlation
    res$covmat_without_rho[,,i_lag] = 
      Matern(h, r = effective_range[i_lag,var_idx], nu_ = nu_[var_idx]) * 
      multiplier[i_lag, var_idx]
    # covariance
    res$covmat[,,i_lag] = res$covmat_without_rho[,,i_lag] * rho_vec_with_ones[var_idx]
  }
  res
}


expand_nu = function(nu_vec){nu_ = outer(nu_vec, nu_vec, "+") ; nu_ = nu_/2; nu_ = nu_[lower.tri(nu_, T)];nu_}

# carry on with toy example

covmat_coeffs_1  =   GMA_compressed(
  locs = locs_1, 
  var_tag = var_tag_1,
  multiplier = multiplier[1,,drop=FALSE], 
  effective_range = effective_range,
  nu_ = expand_nu(nu_vec),
  n_var =  n_var, 
  rho_vec_with_ones = put_ones_in_rho_vec(rho_vec)
)
covmat_coeffs_2  =   GMA_compressed(
  locs = locs_2, 
  var_tag = var_tag_2,
  multiplier = multiplier[-nrow(multiplier),,drop=FALSE], 
  effective_range = effective_range,
  nu_ = expand_nu(nu_vec),
  n_var =  n_var, 
  rho_vec_with_ones = put_ones_in_rho_vec(rho_vec)
)

covmat_coeffs_12  =   GMA_rectangular(
  locs_1 = locs_1, 
  locs_2 = locs_2, 
  var_tag_1 = var_tag_1,
  var_tag_2 = var_tag_2,
  multiplier = multiplier[-1,,drop=FALSE], 
  effective_range = effective_range,
  nu_ = expand_nu(nu_vec),
  n_var =  n_var, 
  rho_vec_with_ones = put_ones_in_rho_vec(rho_vec)
)

# block_lower_tri_idx = lower_tri_idx(n_loc, T)

# expand each sub matrix
expand_covmat_into_blocks = function(covmat_coeffs, n_loc, block_lower_tri_idx)
{
  #if(!sparse)
  blocks = lapply(seq(ncol(covmat_coeffs)), function(i)
  {
    M = matrix(0, n_loc, n_loc)
    M[block_lower_tri_idx]= covmat_coeffs[, i]
    M = t(M)
    M[block_lower_tri_idx]= covmat_coeffs[, i]
    M
  })
  blocks  
}


expand_block_toeplitz_covmat = function(covmat_coeffs, block_lower_tri_idx, n_loc)
{
  blocks = expand_covmat_into_blocks(covmat_coeffs, n_loc, block_lower_tri_idx)
  block_idx_matrix = toeplitz(c(seq(ncol(covmat_coeffs))))
  res=  matrix(0, n_loc*ncol(covmat_coeffs), n_loc*ncol(covmat_coeffs))
  for(i in seq(ncol(block_idx_matrix)))
  {
    for(j in seq(ncol(block_idx_matrix))){
      res[seq(((i-1)*n_loc+1), (i)*n_loc), seq(((j-1)*n_loc+1),(j)*n_loc)]=blocks[[block_idx_matrix[i,j]]]
    }
  }
  res
}

covmat_previous_periods = (expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_2$covmat, n_loc = 50, block_lower_tri_idx = lower_tri_idx(50, T)))
covmat_current_period = expand_covmat_into_blocks(covmat_coeffs = covmat_coeffs_1$covmat, n_loc = 40, block_lower_tri_idx = lower_tri_idx(40, T))[[1]]
side_blocks_rectangles = covmat_coeffs_12$covmat


expand_full_covmat = function(covmat_previous_periods = NULL, covmat_current_period, side_blocks_rectangles = NULL)
{
  if(is.null(covmat_previous_periods)&is.null(side_blocks_rectangles))return(covmat_current_period)
  res = matrix(0, ncol(covmat_previous_periods)+ncol(covmat_current_period), ncol(covmat_previous_periods)+ncol(covmat_current_period))
  res[-seq(ncol(covmat_current_period)), -seq(ncol(covmat_current_period))] = covmat_previous_periods
  res[seq(ncol(covmat_current_period)), seq(ncol(covmat_current_period))] = covmat_current_period
  
  res[seq(ncol(covmat_current_period)), -seq(ncol(covmat_current_period))] = unlist(side_blocks_rectangles)
  res[-seq(ncol(covmat_current_period)), seq(ncol(covmat_current_period))] = t(res[seq(ncol(covmat_current_period)), -seq(ncol(covmat_current_period))])
  res
}





multiply_vector_block_toeplitz_sparse = function(v, covmat_coeffs, n_loc, idx_mat)
{
  k = length(blocks)
  v_ = matrix(v, ncol = k)
  res = matrix(0, n_loc, k)
  for(lag_idx in seq(k))
  {
    toeplitz_idx = rep(0, k)
    toeplitz_idx[lag_idx]=1
    res = res + Multiply_matrix_compressed_symmat(v_, covmat_coeffs[,lag_idx], idx_mat, n_loc) %*% toeplitz(toeplitz_idx)
  }
  res
}

multiply_vector_full_covmat_sparse


#  image(expand_block_toeplitz_covmat(covmat_coeffs$covmat, n_loc, block_lower_tri_idx))
#  try(chol(expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs$covmat, n_loc = n_loc, block_lower_tri_idx = block_lower_tri_idx)))
#  print(diag(expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs$covmat, n_loc = n_loc, block_lower_tri_idx = block_lower_tri_idx)))



