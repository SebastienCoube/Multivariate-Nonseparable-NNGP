


Rcpp::sourceCpp("Scripts/multiply.cpp")
source("Scripts/GMA.R")
source("Scripts/grouping.R")
source("Scripts/multivariate_NN.R")


set.seed(1)

n_loc = 3000
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
                            m_whatever_closest = 50, m_same_var = 5, m_other_vars = 0, 
                            lonlat = F)
DAG = guinness_groups(DAG)
summary(sapply(DAG$children, length))
summary(sapply(DAG$parents, length))



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
  res$
  # looping over each block of the DAG
  for(i in seq(legth(DAG)))
  {
  }
}

