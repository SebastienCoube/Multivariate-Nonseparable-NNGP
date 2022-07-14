gamma0 = function(u,param.gamma0,bb){
  # Returns a vector of length U
  
  # u is a time lag
  # param.gamma0 = (c,a)
  # bb is separability parameter
  cc = param.gamma0[1]
  aa = param.gamma0[2]
  if( (bb < 0) | (aa < 0) | (aa > 1) | (cc<0) ){
    stop("gamma0: parameters out of bounds")} 
  out = (1 + cc*abs(u)^(2*aa))^bb - 1
  return(out)
}

Rcov.uij = function(u,i,j,param.Rcov,A,b){
  # Computes the temporal covariance for one time lag and (i,j)
  
  # u is a time lag
  # param.Rcov = (r,lambda) 
  # b is separability parameter
  # A is the p-vector of weights
  
  r = param.Rcov[1]
  lambda = param.Rcov[2]
  if( (r<0) | (lambda <=0) | (lambda > 1) | (b<0) | (b>1) | any(A<0) | any(A>=1)) stop("Rcov: parameters out of bounds") 
  cv = ( 1 + abs(r*u)^(2*lambda) )^(-b)
  Sij = A[i]*A[j]
  return(cv*Sij)
}

Rcov = function(u,param.Rcov,A,b){
  # Computes the temporal covariance for all time lags and all variables
  # Returns an array [U x p x p]
  
  # u is a time lag
  # param.Rcov = (r,lambda) 
  # b is separability parameter
  # A is the p-vector of weights
  
  r = param.Rcov[1]
  lambda = param.Rcov[2]
  if( (r<0) | (lambda <=0) | (lambda > 1) | (b<0) | any(A<0) | any(A>=1)) stop("Rcov: parameters out of bounds") 
  cv = ( 1 + abs(r*u)^(2*lambda) )^(-b)
  S = A%*%t(A)
  R = array(dim=c(length(u),length(A),length(A)))
  R[1:length(u),,] = cv
  for (tt in 1:length(u)){
    R[tt,1:length(A),1:length(A)] = R[tt,1:length(A),1:length(A)]*S
  }
  return(R)
}

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

MultiGM.hu = function(h,u,rho,param.gamma0,param.Rcov,a,nu,delta,A,b){
  # Computes the full covariance matrix of spatio-temporal function at (h,u)
  # Returns an array [p x p]
  
  p = length(A)
  eps = 10^{-9}
  
  OK = (dim(a)[1]==p) & (dim(a)[2]==p) & (dim(nu)[1]==p) & (dim(nu)[2]==p)
  if (!OK) stop("MulitGM.hu: dimensions don't match")
  if( any(nu < 0) | any(a < 0)  | any(A < 0))  stop("MultiGM.hu: elements of nu, a and A must be postive")
  if ( (b < 0) | (b > 1) ) stop("MultiGM.hu: b should be between 0 and 1")
  if (delta < 0) stop("MultiGM: delta should be positive")
  if (delta > 1) stop("MultiGM: delta should be less or equal to 1")
  
  M = exp(-a)
  if (any(eigen(M)$values < -eps) ) stop("MultiGM.hu: a should be cond. neg. definite")
  M = exp(-nu)
  if (any(eigen(M)$values < -eps) ) stop("MultiGM.hu: nu should be cond. neg. definite")
  M = exp(-a^2+nu)
  if (any(eigen(M)$values < -eps) ) stop("MultiGM.hu: a^2-nu should be cond. neg. definite")
  if (any(A^2>=1))  stop("MultiGM.hu: A_i must less than 1 for all i")
  
  # h is space lag; u is a time lag
  # param.gamma0 = (c,a) as in Eq. (24a) in paper
  # param.Rcov = (r,lambda) as in Eq. (24b) in paper
  # b is separability parameter
  # A is the p-vector of weights
  # a is the p x p range matrix; must be cond. neg. definite
  # nu is the p x p smoothness matrix; must be cond. neg. definite
  # xi is an overall covariance matrix
  # space dimension d = 2 => d/2 = 1
  
  tauvec = diag(a)^{2*diag(nu)}*(1-A^2)^2
  tau = rho*sqrt(outer(tauvec,tauvec,FUN = "*"))
  xi = tau*a^{-2*nu}
  denom =  gamma0(u,param.gamma0,b) + 1 - Rcov(u,param.Rcov,A,b)   # denom without any exponent (since d/2 = 1)
  denom = denom[1,,]
  denom_temp  = gamma0(u,param.gamma0,delta) + 1 - Rcov(u,param.Rcov,A,delta) # denom of parameterization added in Eq. 24
  denom_temp = denom_temp[1,,]
  M = array(dim=c(p,p))
  for (m in 1:p){
    for (n in m:p){
      M[m,n] = xi[m,n]*geoR::matern(h/(denom[m,n])**(1/2),phi=1/a[m,n],kappa = nu[m,n])/
        (denom[m,n]*denom_temp[m,n])
      if(n>m) M[n,m] = M[m,n]
    }
  }
  return(M)
}

MultiGM.huij = function(h,u,i,j,rho,param.gamma0,param.Rcov,a,nu,delta,A,b){
  # Computes the full covariance matrix of spatio-temporal function
  # Returns a single value
  
  p = length(A)
  eps = 10^{-9}
  
  OK = (dim(a)[1]==p) & (dim(a)[2]==p) & (dim(nu)[1]==p) & (dim(nu)[2]==p)
  if (!OK) stop("MulitGM.huij: dimensions don't match")
  if( any(nu < 0) | any(a < 0)  | any(A < 0))  stop("MultiGM.huij: elements of nu, a and A must be postive")
  if ( (b < 0) | (b > 1) ) stop("MultiGM.huij: b should be between 0 and 1")
  if (delta < 0) stop("MultiGM: delta should be positive")
  if (delta > 1) stop("MultiGM: delta should be less or equal to 1")
  
  M = exp(-a)
  if (any(eigen(M)$values < -eps) ) stop("MultiGM.huij: a should be cond. neg. definite")
  M = exp(-nu)
  if (any(eigen(M)$values < -eps) ) stop("MultiGM.huij: nu should be cond. neg. definite")
  M = exp(-a^2+nu)
  if (any(eigen(M)$values < -eps) ) stop("MultiGM.huij: a^2-nu should be cond. neg. definite")
  if (any(A^2>=1))  stop("MultiGM.huij: A_i must less than 1 for all i")
  tauvec = diag(a)^{2*diag(nu)}*(1-A^2)^2
  tau = rho[i,j]*sqrt(tauvec[i]*tauvec[j])
  
  # h is space lag; u is a time lag
  # param.gamma0 = (c,a) 
  # param.Rcov = (r,lambda)
  # b is separability parameter
  # A is the p-vector of weights
  # a is the p x p range matrix; must be cond. neg. definite
  # nu is the p x p smoothness matrix; must be cond. neg. definite
  # xi is an overall covariance matrix
  # space dimension d = 2 => d/2 = 1
  
  xi = tau*a[i,j]^{-2*nu[i,j]}
  denom =  gamma0(u,param.gamma0,b) + 1 - Rcov.uij(u,i,j,param.Rcov,A,b)   # denom without any exponent (since d/2 = 1)
  denom_temp  = gamma0(u,param.gamma0,delta) + 1 - Rcov.uij(u,i,j,param.Rcov,A,delta) # denom of parameterization added in Eq. 24
  X = xi*geoR::matern(h/(denom**(1/2)),phi=1/a[i,j],kappa = nu[i,j])/(denom*denom_temp)
  return(X)
}

MultiGM.0u = function(u,c.hat,param.gamma0,param.Rcov,delta,A,b){
  # Computes the full covariance matrix of temporal function at u
  # Returns an array [p x p]
  
  if ( (b < 0) | (b > 1) ) stop("MultiGM.Ou: b should be between 0 and 1")
  if (delta < 0) stop("MultiGM: delta should be positive")
  if (delta > 1) stop("MultiGM: delta should be less or equal to 1")
  if ( any(A^2>=1) )  stop("MultiGM.ou: A_i must less than 1 for all i")
  
  # u is a time lag
  # c.hat is the experimental p x p covariance matrix
  # param.gamma0 = (c,a)
  # param.Rcov = (r,lambda)
  # b is separability parameter
  # A is the p-vector of weights
  # space dimension d = 2 => d/2 = 1
  
  denom =  gamma0(u,param.gamma0,b) + 1 - Rcov(u,param.Rcov,A,b)   # denom without any exponent (since d/2 = 1)
  denom = denom[1,,]
  denom_temp  = gamma0(u,param.gamma0,delta) + 1 - Rcov(u,param.Rcov,A,delta) # denom of parameterization added in Eq. 24
  denom_temp = denom_temp[1,,]
  M = c.hat*(1-A%*%t(A))^2 / (denom*denom_temp)
  return(M)
}

MultiGM.h0 = function(h,c.hat,a,nu,A){
  # Computes the full covariance matrix of spatial function at h
  # Returns an array [p x p]
  
  p = length(A)
  eps = 10^{-9}
  
  # h is the spatial lag
  # c.hat is the experimental p x p covariance matrix
  # a and nu are p x p matrices
  # A is the p-vector of weights
  
  OK = (dim(a)[1]==p) & (dim(a)[2]==p) & (dim(nu)[1]==p) & (dim(nu)[2]==p)
  if (!OK) stop("MulitGM.h0: dimensions don't match")
  if( any(nu < 0) | any(a < 0)  | any(A < 0))  stop("MultiGM.h0: elements of nu, a and A must be postive")
  M = exp(-a^2+nu)
  if (any(eigen(M)$values < -eps) ) stop("MultiGM.h0: a^2-nu should be cond. neg. definite")
  
  phi = (1-A%*%t(A))^(1/2)/a
  M = array(dim=c(p,p))
  for (m in 1:p){
    for (n in m:p){
      M[m,n] = geoR::matern(h,phi=phi[m,n],kappa = nu[m,n])
      if(n>m) M[n,m] = M[m,n]
    }
  }
  return(c.hat*M)
}

