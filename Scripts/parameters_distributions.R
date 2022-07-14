
beta_non_centered_distribution = function(
    X, # covariates at each site
    z, # observations at each site; can have NAs
    tau, # noise covmats at each sites; dimension may vary ollowing the number of observations
    beta_prior_mean, # prior for the regression coefficients
    beta_prior_precision # prior for the regression coefficients
)
{
  # Amelioration possible : utiliser la formule de Kronecker pour avoir
  # (chol tau z X) T (chol tau z X) =  (chol tau T z XT)  (chol tau z X) = tau z XTX
  
  # for large X, go for Gibbs strategy by dissociating between orthogonal Xes
  
  # passing the noise precisions to Cholesky and transpose
  noise_chols = lapply(tau, function(x)chol(x))
  # weight each observation of X using noise Choleskys
  weighted_X = do.call(rbind,
                       mapply(SIMPLIFY = F,
                              function(row_from_X, row_from_x, noise_precision_chol)
                              {
                                M = matrix(0, length(row_from_x), length(row_from_x))
                                M[!is.na(row_from_x), !is.na(row_from_x)] = noise_precision_chol
                                return(M%x%matrix(row_from_X, nrow = 1))
                              }, 
                              split(X, row(X)), 
                              split(z, row(z)),
                              noise_chols))
  # weight each observation of z using noise Choleskys
  chol_noise_x = do.call(c,
                         mapply(SIMPLIFY = F,
                                function(row_from_x, noise_precision_chol)
                                {
                                  M = matrix(0, length(row_from_x), length(row_from_x))
                                  M[!is.na(row_from_x), !is.na(row_from_x)] = noise_precision_chol
                                  row_from_x_ = row_from_x;row_from_x_[is.na(row_from_x_)] = 0
                                  return(M %*% row_from_x_)
                                }, 
                                split(z, row(z)),
                                noise_chols))
  beta_posterior_precision = beta_prior_precision + crossprod(weighted_X)
  beta_posterior_mean = beta_prior_mean + 
    solve(beta_posterior_precision, crossprod(weighted_X, chol_noise_x))
  
beta_posterior_precision_ = beta_prior_precision 
for(i in seq(nrow(z))){
  print(i)
  M = matrix(0, ncol(z), ncol(z))
  M[!is.na(z[i,]), !is.na(z[i,])] = tau[[i]]
  beta_posterior_precision_ = beta_posterior_precision_ + M %x% tcrossprod(X[i,])
}
  return(list("mean" = beta_posterior_mean, "precision" = beta_posterior_precision))
}