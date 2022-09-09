
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
  t1 = Sys.time()
  # space distance
  space_locs_unique = unique(split(space_locs, row(space_locs)))
  space_locs_match = match(split(space_locs, row(space_locs)), space_locs_unique)
  space_locs_unique = do.call(rbind, space_locs_unique)
  #h = as.matrix(dist(space_locs, diag = T, upper = T))
  h = as.matrix(dist(space_locs_unique, diag = T, upper = T))
  
  # time distance
  time_locs_unique = unique(split(time_locs, row(time_locs)))
  time_locs_match = match(split(time_locs, row(time_locs)), time_locs_unique)
  time_locs_unique = do.call(rbind, time_locs_unique)
  #u = as.matrix(dist(time_locs, diag = T, upper = T))
  u = as.matrix(dist(time_locs_unique, diag = T, upper = T))
  
  
  # pseudo - variogram like in 19, with specification of 5.2
  eta = (
    ((1 + (gamma_0_c*u)^(2*gamma_0_a))^gamma_0_b)[time_locs_match, time_locs_match] # gamma 0 like in 5.2 + 19
    - outer(A[var_tag], A[var_tag]) * ((1 + (R_r * u)^(2*R_lambda))^(-R_b))[time_locs_match, time_locs_match] # Rij like in 5.2 + 19
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
      * (2^(1-nu_mat))[var_tag, var_tag] * (r * h[space_locs_match, space_locs_match])^nu_mat[var_tag, var_tag] * besselK(r * h[space_locs_match, space_locs_match], nu_mat[var_tag, var_tag])# * (1/gamma(nu_mat))[var_tag, var_tag]
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

t1 = Sys.time()
tatat0 = lapply(seq(1000), function(i)
allard_multi_GM(space_locs = space_locs, time_locs = time_locs, var_tag = var_tag, 
                a2 = a2, b2 = b2, nu = nu, alpha = alpha, beta = beta, rho = rho, 
                gamma_0_c = gamma_0_c, gamma_0_a = gamma_0_a, gamma_0_b = gamma_0_b, 
                R_c = R_c, R_lambda = R_lambda, R_b = R_b, A = A
                  )
)
print(Sys.time()-t1)

gamma_0_fun = function(gamma_0_a, gamma_0_b, gamma_0_c, u) (1 + (gamma_0_c*u)^(2*gamma_0_a))^gamma_0_b
 




# first location is conditionned by the other
# slightly modified version of the multivariate Gneiting Matérn covariance of Allard et al 
# https://hal.archives-ouvertes.fr/hal-03564931
# function of equation 19 with some hypotheses from 5.2
# innovation wrt the article: the covariance is scaled on the diagonal so that alpha, beta, 
# a2, b2, nu have no impact on marginal variance

allard_multi_GM_lag = function(
    h, # spatial  lag
    u, # temporal lag
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
  n = length(var_tag)
  
  t1 = Sys.time()
  
  sapply(seq(100), function(x){
    
  # space distance
  space_locs_dist = c(dist(space_locs, diag = T)); space_locs_dist[space_locs_dist==0] = min(space_locs_dist)/1000
  h = unique(space_locs_dist)
  h_to_all = match(space_locs_dist, h)
  # time distance
  time_locs_dist = c(dist(time_locs, diag = T));  time_locs_dist[time_locs_dist==0] = min(time_locs_dist)/1000
  u = unique(time_locs_dist)
  u_to_all = match(time_locs_dist, u)
  # var tag
  var_1 = outer(rep(1, n), var_tag); var_1 = var_1[lower.tri(var_1)]
  var_2 = outer(var_tag, rep(1, n)); var_2 = var_2[lower.tri(var_2)]
  var_tag_combination = paste(pmin(var_1, var_2), pmax(var_1, var_2), sep = "_")
  cmb = unique(var_tag_combination)
  cmb_to_all = match(var_tag_combination, cmb)
  })
  
  print(Sys.time()-t1)
  # crossed indices
  time_vartag = cbind(time_locs_dist, var_tag_combination)
  
  
  
  # time distance
  time_locs_unique = unique(split(time_locs, row(time_locs)))
  time_locs_match = match(split(time_locs, row(time_locs)), time_locs_unique)
  time_locs_unique = do.call(rbind, time_locs_unique)
  #u = as.matrix(dist(time_locs, diag = T, upper = T))
  u = as.matrix(dist(time_locs_unique, diag = T, upper = T))
  u[u==0] = min(u[u!=0])/10000 # dealing with Inf in besselK
  matrix(
  u[
  matrix(match(u, unique(u)), nrow(u))[time_locs_match, time_locs_match]
  ], nrow(time_locs)
  ) - as.matrix(dist(time_locs))
  
  # gamma 0 like in 5.2 + 19
  gamma_0 = (1 + (gamma_0_c*u)^(2*gamma_0_a))^gamma_0_b
  R_0 = outer(A, A) # R^0  like in 5.2 + 19
  R_u_factor = (1 + (R_r * u)^(2*R_lambda))^(-R_b) # Rij / R0 like in 5.2 + 19
  # pseudo - variogram like in 19, with specification of 5.2
  eta = (
    gamma_0[time_locs_match, time_locs_match]
    - R_0[var_tag, var_tag] * (1- R_u_factor[time_locs_match, time_locs_match]) # Rij like in 5.2 + 19
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
  sig_normalization = 
    (
      (1/sqrt(gamma(nu)))* # normalizing diag of gamma(nu). When nu_i \approx nu_j, close to 1. 
      sqrt((alpha + a2)^nu)  * #diag(eta) = 1
      sqrt((beta + b2)^(ncol(space_locs)/2))
    )[var_tag] 
  # test if sig normalization works (* gamma(nu_mat) [var_tag, var_tag] must be uncommented)
  # diag(diag(sig_normalization) %*% sig %*% diag(sig_normalization))
  # matern covariance
  # 1/gamma(nu) is removed because cancels out with previous gamma (nu)
  matern_cov =
    diag(sig_normalization) %*%
    (
      sig 
      * (2^(1-nu_mat))[var_tag, var_tag] * (r * h[space_locs_match, space_locs_match])^nu_mat[var_tag, var_tag] * besselK(r * h[space_locs_match, space_locs_match], nu_mat[var_tag, var_tag])# * (1/gamma(nu_mat))[var_tag, var_tag]
    ) %*%
    diag(sig_normalization)
  # checking that marginal variance corresponds to rho's
  # diag(matern_cov)/diag(rho)[var_tag]
  # checking positive definiteness
  # chol(matern_cov)
  return(matern_cov)
}




allard_multi_GM_vecchia = function(
    space_locs, # spatial locations, first location conditions upon the other
    time_locs,# temporal locations, first location conditions upon the other
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
  space_locs_unique = unique(split(space_locs, row(space_locs)))
  space_locs_match = match(split(space_locs, row(space_locs)), space_locs_unique)
  space_locs_unique = do.call(rbind, space_locs_unique)
  #h = as.matrix(dist(space_locs, diag = T, upper = T))
  h = dist(space_locs_unique)
  h = as.matrix(dist(space_locs_unique, diag = T, upper = T))
  h[h==0] = min(h[h!=0])/10000 # dealing with Inf in besselK
  # time distance
  time_locs_unique = unique(split(time_locs, row(time_locs)))
  time_locs_match = match(split(time_locs, row(time_locs)), time_locs_unique)
  time_locs_unique = do.call(rbind, time_locs_unique)
  #u = as.matrix(dist(time_locs, diag = T, upper = T))
  u = as.matrix(dist(time_locs_unique, diag = T, upper = T))
  u[u==0] = min(u[u!=0])/10000 # dealing with Inf in besselK
  matrix(
  u[
  matrix(match(u, unique(u)), nrow(u))[time_locs_match, time_locs_match]
  ], nrow(time_locs)
  ) - as.matrix(dist(time_locs))
  
  # gamma 0 like in 5.2 + 19
  gamma_0 = (1 + (gamma_0_c*u)^(2*gamma_0_a))^gamma_0_b
  # gamma 0 like in 5.2 + 19
  d_gamma_0_d_c = # derivative
    gamma_0_b * (1 + (gamma_0_c*u)^(2*gamma_0_a))^(gamma_0_b-1) * 
    ((2*gamma_0_a) * (gamma_0_c*u)^(2*gamma_0_a - 1)) *
    u
  # check derivative
  # gamma_0_ = (1 + ((gamma_0_c + 0.0000001)*u)^(2*gamma_0_a))^gamma_0_b; print(1000000 * (gamma_0_ - gamma_0))/d_gamma_0_d_c
  d_gamma_0_d_c2 = # derivative
    (gamma_0_b-1) * gamma_0_b * (1 + (gamma_0_c*u)^(2*gamma_0_a))^(gamma_0_b-2) * 
    ((2*gamma_0_a)^2 * (gamma_0_c*u)^(2 * (2*gamma_0_a - 1))) *
    u^2 + 
    gamma_0_b * (1 + (gamma_0_c*u)^(2*gamma_0_a))^(gamma_0_b-1) * 
    ((2*gamma_0_a) * (2*gamma_0_a - 1) * (gamma_0_c*u)^(2*gamma_0_a - 2)) *
    u^2
  #check derivative
  #print(1000 * (
  #  1000000 * ((1 + ((gamma_0_c + 0.001001)*u)^(2*gamma_0_a))^gamma_0_b - (1 + ((gamma_0_c + 0.001)*u)^(2*gamma_0_a))^gamma_0_b)-
  #    1000000 * ((1 + ((gamma_0_c + 0.000001)*u)^(2*gamma_0_a))^gamma_0_b - (1 + ((gamma_0_c + 0.0)*u)^(2*gamma_0_a))^gamma_0_b)
  #)/d_gamma_0_d_c2)
  t1 = Sys.time()
  tatato = sapply(seq(10000), function(x)
  d_gamma_0_d_c_d_a = # derivative
    (gamma_0_b-1) * gamma_0_b * (1 + (gamma_0_c*u)^(2*gamma_0_a))^(gamma_0_b-2) *
    2 * log((gamma_0_c*u)) * (gamma_0_c*u)^(2*gamma_0_a) *
    ((2*gamma_0_a) * (gamma_0_c*u)^(2*gamma_0_a - 1)) *
    u + 
    gamma_0_b * (1 + (gamma_0_c*u)^(2*gamma_0_a))^(gamma_0_b-1) * 
    (2 * (gamma_0_c*u)^(2*gamma_0_a - 1)) *
    u +
    gamma_0_b * (1 + (gamma_0_c*u)^(2*gamma_0_a))^(gamma_0_b-1) * 
    ((2*gamma_0_a) * 2 * log((gamma_0_c*u)) * (gamma_0_c*u)^(2*gamma_0_a - 1)) *
    u
    )
  print(Sys.time()-t1)
  
  t1 = Sys.time()
  tatato = sapply(seq(10000), function(x)
    tatato = 1000 * (
      1000 * ((1 + ((gamma_0_c + 0.001)*u)^(2*gamma_0_a + 0.001))^gamma_0_b - (1 + ((gamma_0_c + 0.001)*u)^(2*gamma_0_a + 0))^gamma_0_b)-
        1000 * ((1 + ((gamma_0_c + 0    )*u)^(2*gamma_0_a + 0.001))^gamma_0_b - (1 + ((gamma_0_c        )*u)^(2*gamma_0_a + 0))^gamma_0_b)
    )
  )
  print(Sys.time()-t1)
  #check derivative
  print(1000 * (
    1000 * ((1 + ((gamma_0_c + 0.001)*u)^(2*gamma_0_a + 0.001))^gamma_0_b - (1 + ((gamma_0_c + 0.001)*u)^(2*gamma_0_a + 0))^gamma_0_b)-
    1000 * ((1 + ((gamma_0_c + 0    )*u)^(2*gamma_0_a + 0.001))^gamma_0_b - (1 + ((gamma_0_c        )*u)^(2*gamma_0_a + 0))^gamma_0_b)
  )/d_gamma_0_d_c_d_a)
  
  d_gamma_0_d_a = 
    gamma_0_b * (1 + (gamma_0_c*u)^(2*gamma_0_a))^(gamma_0_b-1) * 
    (gamma_0_c*u)^(2*gamma_0_a)

  
  
  
  R_0 = outer(A, A) # R^0  like in 5.2 + 19
  R_u_factor = (1 + (R_r * u)^(2*R_lambda))^(-R_b) # Rij / R0 like in 5.2 + 19
  # pseudo - variogram like in 19, with specification of 5.2
  eta = (
    gamma_0[time_locs_match, time_locs_match]
    - R_0[var_tag, var_tag] * (1- R_u_factor[time_locs_match, time_locs_match]) # Rij like in 5.2 + 19
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
  sig_normalization = 
    (
      (1/sqrt(gamma(nu)))* # normalizing diag of gamma(nu). When nu_i \approx nu_j, close to 1. 
      sqrt((alpha + a2)^nu)  * #diag(eta) = 1
      sqrt((beta + b2)^(ncol(space_locs)/2))
    )[var_tag] 
  # test if sig normalization works (* gamma(nu_mat) [var_tag, var_tag] must be uncommented)
  # diag(diag(sig_normalization) %*% sig %*% diag(sig_normalization))
  # matern covariance
  # 1/gamma(nu) is removed because cancels out with previous gamma (nu)
  matern_cov =
    diag(sig_normalization) %*%
    (
      sig 
      * (2^(1-nu_mat))[var_tag, var_tag] * (r * h[space_locs_match, space_locs_match])^nu_mat[var_tag, var_tag] * besselK(r * h[space_locs_match, space_locs_match], nu_mat[var_tag, var_tag])# * (1/gamma(nu_mat))[var_tag, var_tag]
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




# TRIED TO DO ALL ANALYTIC DIFFERENTIATION
# BUT TOO MUCH FLOPS WRT FINITE DIFF

# slightly modified version of the multivariate Gneiting Matérn covariance of Allard et al 
# https://hal.archives-ouvertes.fr/hal-03564931
# function of equation 19 with some hypotheses from 5.2
# innovation wrt the article: the covariance is scaled on the diagonal so that alpha, beta, 
# a2, b2, nu have no impact on marginal variance

allard_multi_GM_vecchia = function(
  space_locs, # spatial locations, first location conditions upon the other
  time_locs,# temporal locations, first location conditions upon the other
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
  space_locs_unique = unique(split(space_locs, row(space_locs)))
  space_locs_match = match(split(space_locs, row(space_locs)), space_locs_unique)
  space_locs_unique = do.call(rbind, space_locs_unique)
  #h = as.matrix(dist(space_locs, diag = T, upper = T))
  h = as.matrix(dist(space_locs_unique, diag = T, upper = T))
  h[h==0] = min(h[h!=0])/10000 # dealing with Inf in besselK
  
  # time distance
  time_locs_unique = unique(split(time_locs, row(time_locs)))
  time_locs_match = match(split(time_locs, row(time_locs)), time_locs_unique)
  time_locs_unique = do.call(rbind, time_locs_unique)
  #u = as.matrix(dist(time_locs, diag = T, upper = T))
  u = as.matrix(dist(time_locs_unique, diag = T, upper = T))
  u[u==0] = min(u[u!=0])/10000 # dealing with Inf in besselK
  
  
  # gamma 0 like in 5.2 + 19
  gamma_0 = (1 + (gamma_0_c*u)^(2*gamma_0_a))^gamma_0_b
  d_gamma_0_d_c = # derivative
    gamma_0_b * (1 + (gamma_0_c*u)^(2*gamma_0_a))^(gamma_0_b-1) * 
    ((2*gamma_0_a) * (gamma_0_c*u)^(2*gamma_0_a - 1)) *
    u
  # check derivative
  # gamma_0_ = (1 + ((gamma_0_c + 0.0000001)*u)^(2*gamma_0_a))^gamma_0_b; print(1000000 * (gamma_0_ - gamma_0))/d_gamma_0_d_c
  d_gamma_0_d_c2 = # derivative
    (gamma_0_b-1) * gamma_0_b * (1 + (gamma_0_c*u)^(2*gamma_0_a))^(gamma_0_b-2) * 
    ((2*gamma_0_a)^2 * (gamma_0_c*u)^(2 * (2*gamma_0_a - 1))) *
    u^2 + 
    gamma_0_b * (1 + (gamma_0_c*u)^(2*gamma_0_a))^(gamma_0_b-1) * 
    ((2*gamma_0_a) * (2*gamma_0_a - 1) * (gamma_0_c*u)^(2*gamma_0_a - 2)) *
    u^2
  #check derivative
  print(1000 * (
    1000000 * ((1 + ((gamma_0_c + 0.001001)*u)^(2*gamma_0_a))^gamma_0_b - (1 + ((gamma_0_c + 0.001)*u)^(2*gamma_0_a))^gamma_0_b)-
    1000000 * ((1 + ((gamma_0_c + 0.000001)*u)^(2*gamma_0_a))^gamma_0_b - (1 + ((gamma_0_c + 0.0)*u)^(2*gamma_0_a))^gamma_0_b)
  )/d_gamma_0_d_c2)
  d_gamma_0_d_c_d_a = # derivative
    gamma_0_b * (1 + (gamma_0_c*u)^(2*gamma_0_a))^(gamma_0_b-1) * 
    ((2*gamma_0_a) * (gamma_0_c*u)^(2*gamma_0_a - 1)) *
    u + 
    gamma_0_b * (1 + (gamma_0_c*u)^(2*gamma_0_a))^(gamma_0_b-1) * 
    ((2*gamma_0_a) * (gamma_0_c*u)^(2*gamma_0_a - 1)) *
    u
  
  R_0 = outer(A, A) # R^0  like in 5.2 + 19
  R_u_factor = (1 + (R_r * u)^(2*R_lambda))^(-R_b) # Rij / R0 like in 5.2 + 19
  # pseudo - variogram like in 19, with specification of 5.2
  eta = (
    gamma_0[time_locs_match, time_locs_match]
    - R_0[var_tag, var_tag] * (1- R_u_factor[time_locs_match, time_locs_match]) # Rij like in 5.2 + 19
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
  # 1/gamma(nu) is removed because cancels out with previous gamma (nu)
  matern_cov =
    diag(sig_normalization) %*%
    (
      sig 
      * (2^(1-nu_mat))[var_tag, var_tag] * (r * h[space_locs_match, space_locs_match])^nu_mat[var_tag, var_tag] * besselK(r * h[space_locs_match, space_locs_match], nu_mat[var_tag, var_tag])# * (1/gamma(nu_mat))[var_tag, var_tag]
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

