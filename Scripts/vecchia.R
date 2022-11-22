


Rcpp::sourceCpp("Scripts/multiply.cpp")
source("Scripts/GMA.R")
source("Scripts/grouping.R")
source("Scripts/multivariate_NN.R")


set.seed(1)

n_loc = 1000
n_var = 5

# spatial locations and variable index
locs = cbind(runif(n_loc), runif(n_loc))
var_tag = ceiling(n_var*runif(n_loc));var_tag = match(var_tag, unique(var_tag))

# covariance params
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

# creating the DAG
DAG = find_ordered_nn_multi(locs = locs, 
                            var_tag = var_tag, 
                            m_whatever_closest = 30, m_same_var = 0, m_other_vars = 0, 
                            lonlat = F)
DAG = guinness_groups(DAG)
summary(sapply(DAG$children, length))
hist(sapply(DAG$children, length), breaks = seq(1, max(sapply(DAG$children, length))))
summary(sapply(DAG$parents, length))


# Pre-computing elements
multiplier = get_multiplier(
  a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
  r = r, A_vec = A_vec, nu_vec = nu_vec, a2_vec = a2_vec, u = u
)
effective_range = get_effective_range(
  a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
  r = r, A_vec = A_vec, a2_vec = a2_vec, u = u
)





i=1

multi_Vecchia = function(
  DAG, 
  locs, var_tag, 
  multiplier, effective_range, 
  rho_vec, 
  compute_derivative_wrt_rho = F
  )
{
  res = list()
  coeffs = list()
  
  n_time = nrow(multiplier)
  n_var = max(var_tag)
  nu_ = expand_nu(nu_vec)
  rho_vec_with_ones = put_ones_in_rho_vec(rho_vec)
  # looping over each block of the DAG
  t1 = Sys.time()
  for(i in seq(length(DAG[[1]])))
  {
    locs_ = locs[c(DAG$children[[i]], DAG$parents[[i]]),]
    var_tag_ = var_tag[c(DAG$children[[i]], DAG$parents[[i]])]
    n_loc_ = nrow(locs_)
    block_lower_tri_idx = lower_tri_idx(n_loc_, diag = T)

    covmat_coeffs = GMA_compressed(
      locs = locs_, 
      var_tag = var_tag_,
      multiplier = multiplier, 
      effective_range = effective_range,
      nu_ = nu_, 
      rho_vec_with_ones = rho_vec_with_ones, 
      n_var = n_var
    )$covmat
    
    covmat_chol = chol(
      expand_block_toeplitz_covmat(
        covmat_coeffs = covmat_coeffs, 
        block_lower_tri_idx = block_lower_tri_idx, 
        n_loc = n_loc_
      )
    )
    
    coeffs[[i]] = 
    backsolve(r =  covmat_chol
      , diag(1, n_loc_*n_time, length(DAG$children[[i]]))
    )
  }
  Sys.time()-t1
  
  
}


