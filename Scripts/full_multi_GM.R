
# slightly modified version of the multivariate Gneiting Matérn covariance of Allard et al 
# https://hal.archives-ouvertes.fr/hal-03564931
# function of equation 19 with some hypotheses from 5.2
# innovation wrt the article: the covariance is scaled on the diagonal so that alpha, beta, 
# a2, b2, nu have no impact on marginal variance

allard_multi_GM = function(
  space_locs, 
  time_locs, 
  var_tag, # variable tag, must be integer
  a2, b2, nu, # positive vectors of size p to be converted in a-separable matrices
  # a2/b2 = spatial range at time 0
  # nu = smoothness
  alpha, beta, # nonnegative numbers
  # alpha/beta = spatial range at infinite time
  rho, # between-variables variance (without the effect of nu) parametrized like in 5.2
  gamma_0_c, # number greater than 0, controlling time correlation
  gamma_0_a, gamma_0_b, # numbers between 0 and 1, controlling time correlation
  R_c, R_lambda, R_b, # numbers between 0 and 1, controlling time correlation
  A # vector of size p with elements between 0 and 1, controlling time correlation
){
  # space distance
  h = as.matrix(dist(space_locs, diag = T, upper = T))
  # time distance
  u = as.matrix(dist(time_locs, diag = T, upper = T))
  # pseudo - variogram like in 19, with specification of 5.2
  eta = (
    (1 + (gamma_0_c*u)^(2*gamma_0_a))^gamma_0_b # gamma 0 like in 5.2 + 19
    - outer(A[var_tag], A[var_tag]) * (1 + (R_r * u)^(2*R_lambda))^(-R_b) # Rij like in 5.2 + 19
    + outer(A[var_tag], A[var_tag]) # R^0  like in 5.2 + 19
  )
  # r like in 15 
  a2_mat = (outer(a2, a2, "+")/2); b2_mat = (outer(b2, b2, "+")/2)
  r = sqrt(
    (alpha * eta + a2_mat[var_tag, var_tag]) / 
    (beta  * eta + b2_mat[var_tag, var_tag])
    )
  # plot(r, u, col = outer(var_tag, var_tag, function(x, y)as.factor(paste(x, y)) ))
  # sigma like in 15
   # removing exp(nu) since exp(-nu) is p-d
   # removing gamma(nu) which cancels out with Matern 1/gamma(nu)
  nu_mat = outer(nu, nu, "+")/2
  sig = rho[var_tag, var_tag] / 
    (
     (alpha * eta + a2_mat[var_tag, var_tag])^nu_mat[var_tag, var_tag]
     *(beta * eta + b2_mat[var_tag, var_tag])^(ncol(space_locs)/2)
    ) # * gamma(nu_mat) [var_tag, var_tag]
  sig_normalization = (1/sqrt(gamma(nu)))[var_tag] * # normalizing diag of gamma(nu). When nu_i \approx nu_j, close to 1. 
    sqrt((alpha + a2)^nu)[var_tag]  * #diag(eta) = 1
    sqrt((beta + b2)^(ncol(space_locs)/2))[var_tag]  #diag(eta) = 1
  # test if sig normalization works (* gamma(nu_mat) [var_tag, var_tag] must be uncommented)
  # diag(diag(sig_normalization) %*% sig %*% diag(sig_normalization))
  # matern covariance
  h[h==0] = min(h[h!=0])/10000 # dealing with Inf in besselK
  # 1/gamma(nu) is removed because cancels out with previous gamma (nu)
  matern_cov =
    diag(sig_normalization) %*%
    (
      sig 
      * (2^(1-nu_mat))[var_tag, var_tag] * (r * h)^nu_mat[var_tag, var_tag] * besselK(r * h, nu_mat[var_tag, var_tag])# * (1/gamma(nu_mat))[var_tag, var_tag]
    ) %*%
    diag(sig_normalization)
  # checking that marginal variance corresponds to rho's
  # diag(matern_cov)/diag(rho)[var_tag]
  # checking positive definiteness
  # chol(matern_cov)
  return(matern_cov)
}

#space_locs = cbind(runif(100), runif(100))
space_locs = cbind(rep(runif(10), 10), rep(runif(10), 10))
time_locs = cbind(rep(seq(10), each = 10))
p = 3
var_tag = rep(1+floor(p*runif(10)), 10)
a2 = 100 * runif(p)
b2 = 100 * runif(p)
nu = .25 + runif(p)
alpha = runif(1)
beta = runif(1)
rho = crossprod(matrix(rnorm(p^2), p))
gamma_0_a = runif(1) 
gamma_0_b = runif(1) 
gamma_0_c = 3*runif(1) 
R_lambda = runif(1) 
R_b = runif(1) 
R_r = 3*runif(1) 
A = runif(p)





allard_multi_GM = function(
    space_locs, 
    time_locs, 
    var_tag, # variable tag, must be integer
    a2, b2, nu, # positive vectors of size p to be converted in a-separable matrices
    # a2/b2 = spatial range at time 0
    # nu = smoothness
    alpha, beta, # nonnegative numbers
    # alpha/beta = spatial range at infinite time
    rho, # between-variables variance (without the effect of nu) parametrized like in 5.2
    gamma_0_c, # number greater than 0, controlling time correlation
    gamma_0_a, gamma_0_b, # numbers between 0 and 1, controlling time correlation
    R_c, R_lambda, R_b, # numbers between 0 and 1, controlling time correlation
    A # vector of size p with elements between 0 and 1, controlling time correlation
){
  # space distance
  h = as.matrix(dist(space_locs, diag = T, upper = T))
  # time distance
  u = as.matrix(dist(time_locs, diag = T, upper = T))
  # pseudo - variogram like in 19, with specification of 5.2
  eta = (
    (1 + (gamma_0_c*u)^(2*gamma_0_a))^gamma_0_b # gamma 0 like in 5.2 + 19
    - outer(A[var_tag], A[var_tag]) * (1 + (R_r * u)^(2*R_lambda))^(-R_b) # Rij like in 5.2 + 19
    + outer(A[var_tag], A[var_tag]) # R^0  like in 5.2 + 19
  )
  # r like in 15 
  a2_mat = (outer(a2, a2, "+")/2); b2_mat = (outer(b2, b2, "+")/2)
  r = sqrt(
    (alpha * eta + a2_mat[var_tag, var_tag]) / 
      (beta  * eta + b2_mat[var_tag, var_tag])
  )
  # plot(r, u, col = outer(var_tag, var_tag, function(x, y)as.factor(paste(x, y)) ))
  # sigma like in 15
  # removing exp(nu) since exp(-nu) is p-d
  # removing gamma(nu) which cancels out with Matern 1/gamma(nu)
  nu_mat = outer(nu, nu, "+")/2
  sig = rho[var_tag, var_tag] / 
    (
      (alpha * eta + a2_mat[var_tag, var_tag])^nu_mat[var_tag, var_tag]
      *(beta * eta + b2_mat[var_tag, var_tag])^(ncol(space_locs)/2)
    ) # * gamma(nu_mat) [var_tag, var_tag]
  sig_normalization = (1/sqrt(gamma(nu)))[var_tag] * # normalizing diag of gamma(nu). When nu_i \approx nu_j, close to 1. 
    sqrt((alpha + a2)^nu)[var_tag]  * #diag(eta) = 1
    sqrt((beta + b2)^(ncol(space_locs)/2))[var_tag]  #diag(eta) = 1
  # test if sig normalization works (* gamma(nu_mat) [var_tag, var_tag] must be uncommented)
  # diag(diag(sig_normalization) %*% sig %*% diag(sig_normalization))
  # matern covariance
  h[h==0] = min(h[h!=0])/10000 # dealing with Inf in besselK
  # 1/gamma(nu) is removed because cancels out with previous gamma (nu)
  matern_cov =
    diag(sig_normalization) %*%
    (
      sig 
      * (2^(1-nu_mat))[var_tag, var_tag] * (r * h)^nu_mat[var_tag, var_tag] * besselK(r * h, nu_mat[var_tag, var_tag])# * (1/gamma(nu_mat))[var_tag, var_tag]
    ) %*%
    diag(sig_normalization)
  # checking that marginal variance corresponds to rho's
  # diag(matern_cov)/diag(rho)[var_tag]
  # checking positive definiteness
  # chol(matern_cov)
  return(matern_cov)
}





multivariate_GM_first = function(
  space_locs, 
  time_locs, 
  var_tag, 
  r2, nu, # positive vectors of size p to be converted in a-separable matrices
  rho, # between-variables correlation (without the effect of nu)
  gamma_0_c, gamma_0_a, gamma_0_b, # numbers between 0 and 1
  R_c, R_lambda, R_b, # numbers between 0 and 1
  A # vector of size p with elements between 0 and 1 
){
  h = as.matrix(dist(space_locs, diag = T, upper = T))
  u = as.matrix(dist(time_locs, diag = T, upper = T))
  nu_mat = outer(nu, nu, "+")/2
  r_mat = sqrt(outer(r2, r2, "+")/2)
  penalty_mat = gamma(nu_mat)/outer(sqrt(gamma(nu)), sqrt(gamma(nu))) * 
    
  
  (
  (1+(gamma_0_c*u)^(2*gamma_0_a))^gamma_0_b - 
  outer(A[var_tag], A[var_tag]) * (1 + (R_r * u)^(2*R_lambda))^(-R_b)
  )
}



space_locs = cbind(runif(100), runif(100))
time_locs = cbind(seq(100))
p = 3
var_tag = 1+floor(p*runif(100));var_tag = match(var_tag, unique(var_tag))
r2 = runif(p)
nu = .25 + runif(p)
alpha = 1
beta = 1
rho = cor(matrix(rnorm(p^2), p))
gamma_0_a = runif(1) 
gamma_0_b = runif(1) 
gamma_0_c = 3*runif(1) 
R_lambda = runif(1) 
R_b = runif(1) 
R_r = 3*runif(1) 
A = runif(p)
