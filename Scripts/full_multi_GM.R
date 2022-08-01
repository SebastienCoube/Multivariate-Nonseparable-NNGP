multivariate_GM = function(
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
   # multiplying by 1/sqrt(diag(gamma(nu))) in order to get unit correlation
  nu_mat = outer(nu, nu, "+")
  #   gamma_nu_normaliztion =  outer(1/sqrt(gamma(2*nu)), 1/sqrt(gamma(2*nu)), "*")
  #   * gamma(nu_mat)[var_tag, var_tag] * gamma_nu_normaliztion[var_tag, var_tag] 
  # removing gamma(nu) because it cancels out with gamma in Matern function
  sig = rho[var_tag, var_tag]  / 
    (
     (alpha * eta + a2_mat[var_tag, var_tag])^nu_mat[var_tag, var_tag]
     *(beta * eta + b2_mat[var_tag, var_tag])^(ncol(space_locs)/2)
    )
  # matern covariance
  # nu being cond n-d, 2^-nu is 
  diag(h) = .0001
  matern_cov =
  sig * (r * h)^nu_mat[var_tag, var_tag] * besselK(r * h, nu_mat[var_tag, var_tag])
  diag( (r * h)^nu_mat[var_tag, var_tag] * besselK(r * h, nu_mat[var_tag, var_tag]))
  diag(sig*gamma(nu_mat[var_tag, var_tag])*exp(nu_mat[var_tag, var_tag]))
  diag(matern_cov)
  chol(matern_cov)
}

space_locs = cbind(runif(100), runif(100))
time_locs = cbind(rep(seq(10), each = 10))
p = 3
var_tag = 1+floor(p*runif(100));var_tag = match(var_tag, unique(var_tag))
a2 = 100 * runif(p)
b2 = 100 * runif(p)
nu = .25 + runif(p)
alpha = runif(1)
beta = runif(1)
rho = diag(1, p)
gamma_0_a = runif(1) 
gamma_0_b = runif(1) 
gamma_0_c = 3*runif(1) 
R_lambda = runif(1) 
R_b = runif(1) 
R_r = 3*runif(1) 
A = runif(p)



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
