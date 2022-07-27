
MultiGM = function(h,u,rho,param.gamma0,param.Rcov,a,nu,delta,A,b){
  # Computes the full covariance matrix of spatio-temporal function
  # Returns an array [H x U x p x p]
  
  p = length(A)
  eps = 10^{-9}
  
  OK = (dim(a)[1]==p) & (dim(a)[2]==p) & (dim(nu)[1]==p) & (dim(nu)[2]==p)
  if (!OK) stop("MulitGM: dimensions don't match")
  if( any(nu < 0) | any(a < 0)  | any(A < 0))  stop("MultiGM: elements of nu, a and A must be postive")
  if ( (b < 0) | (b > 1) ) stop("MultiGM: b should be between 0 and 1")
  if (delta < 0) stop("MultiGM: delta should be positive")
  M = exp(-a)
  if (any(eigen(M)$values < -eps) ) stop("MultiGM: a should be cond. neg. definite")
  M = exp(-nu)
  if (any(eigen(M)$values < -eps) ) stop("MultiGM: nu should be cond. neg. definite")
  M = exp(-a^2+nu)
  if (any(eigen(M)$values < -eps) ) stop("MultiGM: a^2-nu should be cond. neg. definite")
  if (any(A^2>=1))  stop("MultiGM: A_i must less than 1 for all i")
  tauvec = diag(a)^{2*diag(nu)}*(1-A^2)^2
  tau = rho*sqrt(outer(tauvec,tauvec,FUN = "*"))
  M = tau*exp(-nu)/gamma(nu)
  if (any(eigen(M)$values < -eps) ) stop("MultiGM: tau*e^{-nu}/Gamma(nu) should be pos. semidefinite (corollary 1)")
  
  # h is space lag; u is a time lag
  # param.gamma0 = (c,a) 
  # param.Rcov = (r,lambda) 
  # b is separability parameter
  # A is the p-vector of weights
  # a is the p x p range matrix; must be cond. neg. definite
  # nu is the p x p smoothness matrix; must be cond. neg. definite
  # xi is an overall covariance matrix
  # space dimension d = 2 => d/2 = 1
  
  xi = tau*a^{-2*nu}
  denom =  gamma0(u,param.gamma0,b) + 1 - Rcov(u,param.Rcov,A,b)   # denom without any exponent (since d/2 = 1)
  denom_temp  = gamma0(u,param.gamma0,delta) + 1 - Rcov(u,param.Rcov,A,delta) # denom of parameterization added in Eq. 24
  M = array(dim=c(length(h),length(u),p,p))
  
  for (s in 1:length(h)){
    for (m in 1:p){
      for (n in m:p){
        M[s,,m,n] = xi[m,n]*geoR::matern(h[s]/(denom[,m,n])**(1/2),phi=1/a[m,n],kappa = nu[m,n])/(denom[,m,n]*denom_temp[,m,n])
        if(n>m) M[s,,n,m] = M[s,,m,n]
      }
    }
  }
  return(M)
}



h = dist(cbind(1, rnorm(40)))
u = dist(cbind(1, rnorm(40)))
rho = matrix(c(1, .5, .5, 1), 2)
a = sqrt(outer(c(1, .5), c(1, .5), "+"))
nu = outer(c(.8, .5), c(.8, .5), "+")
b = .5
A = c(.8, .5)
delta = 1
param.gamma0 = c(.5, .5)
param.Rcov = c(2, .5)
