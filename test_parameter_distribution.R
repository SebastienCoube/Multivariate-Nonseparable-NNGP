
source("parameters_distributions.R")

X = matrix(rnorm(10000000), 10000)
beta = matrix(rnorm(2000), 1000)
z = X %*% beta
#z[seq(100), 1] = NA
#z[seq(101, 200), 2] = NA
n_obs = ncol(z)-apply(z, 1, function(x)sum(is.na(x)))
tau = lapply(n_obs, function(x)
{
  tatato = rnorm(x*(x+1)/2, 0)
  M = matrix(0,x,x)
  M[lower.tri(M, diag = T)] = tatato
  M = t(M)
  M[lower.tri(M, diag = T)] = tatato
  expm::expm(M)
}
)
z = t(z)
z[!is.na(z)] = z[!is.na(z)] +
  unlist(lapply(tau, function(tau)t(chol(tau))%*%rnorm(nrow(tau))))
z = t(z)

beta_prior_precision = diag(0, length(beta), length(beta))
beta_prior_mean = rep(0, length(beta))

test = beta_non_centered_distribution (
  X = X, # covariates at each site
  z = z, # observations at each site; can have NAs
  tau = tau, # noise covmats at each sites; dimension may vary ollowing the number of observations
  beta_prior_mean = beta_prior_mean, # prior for the regression coefficients
  beta_prior_precision  = beta_prior_precision # prior for the regression coefficients
)
crossprod(X, rnorm(10000)*X)
