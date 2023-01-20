#################
# Generate data #
#################

#### # flip image function to match matrix
#### my_image = function(m)image(t(m)[,nrow(m):1])
#### 
#### # generate test data
#### set.seed(10)
#### n_var = 2
#### 
#### locs_1 = 10*cbind(runif(40), runif(40))
#### locs_2 = rbind(locs_1, 10*cbind(runif(50), runif(50)))
#### var_tag_1 = 1+floor(n_var*runif(40))
#### var_tag_2 = c(var_tag_1, 1+floor(n_var*runif(50)))
#### 
#### rho = GpGp::exponential_isotropic(c(1, 1, 0), .1*matrix(rnorm(2*n_var), n_var))
#### rho_vec = rho[lower.tri(rho, diag = F)];  remove(rho)
#### 
#### a2_vec = 1*runif(n_var)
#### nu_vec = .5 + 2*runif(n_var)
#### alpha = 1 * runif(1)
#### a = runif(1) 
#### b = runif(1) 
#### cc = .1 * runif(1) 
#### lambda = runif(1) 
#### delta = runif(1) 
#### r = .1 * runif(1) 
#### A_vec = runif(n_var)
#### u = seq(0, 5)


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

#### diag_indices_lower_tri = function(n)cumsum(c(1,seq(n, 2)))
#### diag_indices_lower_tri(3)
#### diag_indices_lower_tri(4)
#### diag_indices_lower_tri(5)

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

get_lower_tri_idx = function(n_loc, diag = F)
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
#### get_lower_tri_idx(10)

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
####  i = 1;j =3;n_loc=3;position_in_lower_tri(i, j, n_loc)
####  i = 2;j =1;n_loc=3;position_in_lower_tri(i, j, n_loc)
####  i = 2;j =3;n_loc=3;position_in_lower_tri(i, j, n_loc)
####  i = 3;j =3;n_loc=3;position_in_lower_tri(i, j, n_loc)


# given a vector of integers i_vec and an integer n_var, 
# n_loc var being greater than the entries of i_vec, 
# computes the lower-triangular without diagonal of the crossing matrix of i_vec, i_vec
# and finds the position of those pairs in the lower triangular part of the 
# matrix of size n_var * n_var
position_in_lower_tri_cross_vec = function(i_vec, n_var, diag = F)
{  
  idx = get_lower_tri_idx(length(i_vec), diag = diag)
  position_in_lower_tri(i_vec[idx[,1]], i_vec[idx[,2]], n_var)
}

####  i_vec = 1 + floor(n_var*runif(10))
####  position_in_lower_tri_cross_vec(i_vec, n_var)
####  length(position_in_lower_tri_cross_vec(i_vec, n_var))
####  
####  i_vec = 1 + floor(n_var*runif(10))
####  position_in_lower_tri_cross_vec(i_vec, n_var)
####  length(position_in_lower_tri_cross_vec(i_vec,n_var))

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

#### # Toy example  
#### multiplier = get_multiplier(
####   a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
####   r = r, A_vec = A_vec, nu_vec = nu_vec, a2_vec = a2_vec, u = u
#### )
#### effective_range = get_effective_range(
####   a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
####   r = r, A_vec = A_vec, a2_vec = a2_vec, u = u
#### )

# compute unique elements of GMA covmat
GMA_compressed = function(
    locs, 
    var_idx,
    lower_tri_idx,
    multiplier, 
    effective_range,
    nu_, 
    rho_vec_with_ones, 
    return_with_rho = F
)
{
  # size of spatial sets of interest
  n_loc = nrow(locs)
  n_lags = nrow(multiplier)
  # distance between locs pairs
  h = as.matrix(dist(locs, diag = T));h = h[lower.tri(h, diag = T)]
  if(all(h==0)) h = .00001
  h[h==0] = min(h[h!=0])*.0001
  # creating matrices
  res = list()
  if(return_with_rho)res$covmat_without_rho = matrix(0, nrow(lower_tri_idx), n_lags)
  res$covmat = matrix(0, nrow(lower_tri_idx), n_lags)
  # index of position in the covariance matrix
  if(return_with_rho)
  {
    # looping over time lags
    for(i_lag in seq(n_lags)){
      # covariance without correlation
      res$covmat_without_rho[,i_lag] = 
        Matern(h, r = effective_range[i_lag,var_idx], nu_ = nu_[var_idx]) * 
        multiplier[i_lag, var_idx]
      # covariance
      res$covmat[,i_lag] = res$covmat_without_rho[,i_lag] * rho_vec_with_ones[var_idx]
    }
  }
  if(!return_with_rho)
  {
    # looping over time lags
    for(i_lag in seq(n_lags)){
      # covariance
      res$covmat[,i_lag] = Matern(h, r = effective_range[i_lag,var_idx], nu_ = nu_[var_idx]) * 
        multiplier[i_lag, var_idx] * rho_vec_with_ones[var_idx]
    }
  }
  res
}

# elements of rectangular GMA matrix
GMA_rectangular = function(
    locs_1, locs_2, 
    var_idx,
    multiplier, 
    effective_range,
    nu_, 
    rho_vec_with_ones, 
    n_var, 
    return_with_rho = F
)
{
  # size of spatial sets of interest
  n_loc_1 = nrow(locs_1)
  n_loc_2 = nrow(locs_1)
  n_lags = nrow(multiplier)
  # distance between locs pairs
  h = fields::rdist(locs_1, locs_2)
  h[h==0] = min(h[h!=0])*.0001
  # creating matrices
  res = list()
  if(return_with_rho) res$covmat_without_rho = array(0, c(dim(h), nrow(multiplier)))
  res$covmat = array(0, c(dim(h), nrow(multiplier)))
  
  if(return_with_rho){
    # looping over time lags
    for(i_lag in seq(n_lags)){
      # index of position in the covariance matrix
      # covariance without correlation
      res$covmat_without_rho[,,n_lags-i_lag+1] = 
        Matern(h, r = effective_range[i_lag,var_idx], nu_ = nu_[var_idx]) * 
        multiplier[i_lag, var_idx]
      # covariance
      res$covmat[,,n_lags-i_lag+1] = res$covmat_without_rho[,,n_lags-i_lag+1] * rho_vec_with_ones[var_idx]
    }
  }
  if(!return_with_rho){
    # looping over time lags
    for(i_lag in seq(n_lags)){
      # covariance
      res$covmat[,,n_lags-i_lag+1] = Matern(h, r = effective_range[i_lag,var_idx], nu_ = nu_[var_idx]) * 
        multiplier[i_lag, var_idx] * rho_vec_with_ones[var_idx]
    }
  }
  res
}


expand_nu = function(nu_vec){nu_ = outer(nu_vec, nu_vec, "+") ; nu_ = nu_/2; nu_ = nu_[lower.tri(nu_, T)];nu_}

#### # carry on with toy example
#### covmat_coeffs_1  =   GMA_compressed(
####   locs = locs_1, 
####   lower_tri_idx = get_lower_tri_idx(nrow(locs_1), diag = T),
####   var_idx = position_in_lower_tri_cross_vec(var_tag_1, n_var, diag = T),
####   multiplier = multiplier[1,,drop=FALSE], 
####   effective_range = effective_range,
####   nu_ = expand_nu(nu_vec),
####   rho_vec_with_ones = put_ones_in_rho_vec(rho_vec)
#### )
#### covmat_coeffs_2  =   GMA_compressed(
####   locs = locs_2, 
####   lower_tri_idx = get_lower_tri_idx(nrow(locs_2), diag = T),
####   var_idx = position_in_lower_tri_cross_vec(var_tag_2, n_var, diag = T),
####   multiplier = multiplier[-nrow(multiplier),,drop=FALSE], 
####   effective_range = effective_range,
####   nu_ = expand_nu(nu_vec),
####   rho_vec_with_ones = put_ones_in_rho_vec(rho_vec)
#### )
#### 
#### 
#### tatato = expand.grid(var_tag_1, var_tag_2)
#### var_idx_12 = position_in_lower_tri(tatato[,1], tatato[,2], n_var)
#### covmat_coeffs_12  =   GMA_rectangular(
####   locs_1 = locs_1, 
####   locs_2 = locs_2, 
####   var_idx = var_idx_12,
####   multiplier = multiplier[-1,,drop=FALSE], 
####   effective_range = effective_range[-1,,drop=FALSE],
####   nu_ = expand_nu(nu_vec),
####   n_var =  n_var, 
####   rho_vec_with_ones = put_ones_in_rho_vec(rho_vec)
#### )

# block_lower_tri_idx = get_lower_tri_idx(n_loc, T)

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

#### covmat_previous_periods = (expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_2$covmat, n_loc = 90, block_lower_tri_idx = get_lower_tri_idx(90, T)))
#### covmat_current_period = expand_covmat_into_blocks(covmat_coeffs = covmat_coeffs_1$covmat, n_loc = 40, block_lower_tri_idx = get_lower_tri_idx(40, T))[[1]]
#### side_blocks_rectangles = covmat_coeffs_12$covmat


expand_full_covmat = function(covmat_previous_periods = NULL, covmat_current_period, side_blocks_rectangles = NULL)
{
  if(is.null(covmat_previous_periods)&is.null(side_blocks_rectangles))return(covmat_current_period)
  res = matrix(0, ncol(covmat_previous_periods)+ncol(covmat_current_period), ncol(covmat_previous_periods)+ncol(covmat_current_period))
  res[seq(ncol(covmat_previous_periods)), seq(ncol(covmat_previous_periods))] = covmat_previous_periods
  res[-seq(ncol(covmat_previous_periods)), -seq(ncol(covmat_previous_periods))] = covmat_current_period
  
  res[-seq(ncol(covmat_previous_periods)), seq(ncol(covmat_previous_periods))] = unlist(side_blocks_rectangles)
  res[seq(ncol(covmat_previous_periods)), -seq(ncol(covmat_previous_periods))] = t(res[-seq(ncol(covmat_previous_periods)), seq(ncol(covmat_previous_periods))])
  res
}

#### full_covmat = expand_full_covmat(covmat_previous_periods = expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_2$covmat, block_lower_tri_idx = get_lower_tri_idx(90, T), n_loc = 90), 
####                                  covmat_current_period   = expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_1$covmat, block_lower_tri_idx = get_lower_tri_idx(40, T), n_loc = 40), 
####                                  side_blocks_rectangles  = covmat_coeffs_12$covmat
####                                  )
#### chol(full_covmat)
#### my_image(log(.5+full_covmat))


