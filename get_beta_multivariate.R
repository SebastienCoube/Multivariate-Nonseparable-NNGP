

X = cbind(rep(1, 1000), seq(1000)*.01 )
beta = matrix(rnorm(6), 2)
x = X %*% beta
x[seq(100), 1] = NA
x[seq(101, 200), 2] = NA
x[seq(201, 300), 3] = NA
x[seq(301, 400), c(1, 2)] = NA
n_obs = ncol(x)-apply(x, 1, function(x)sum(is.na(x)))
tau = lapply(n_obs, function(x)
{
  tatato = rnorm(x*(x+1)/2, +1)
  M = matrix(0,x,x)
  M[lower.tri(M, diag = T)] = tatato
  M = t(M)
  M[lower.tri(M, diag = T)] = tatato
  expm::expm(M)
}
)
x = t(x)
x[!is.na(x)] = x[!is.na(x)] +
unlist(lapply(tau, function(tau)t(chol(solve(tau)))%*%rnorm(nrow(tau))))
x = t(x)

beta_prior_precision = diag(.000001, 6,6)
beta_prior_mean = c(0, 0, 0, 0, 0, 0)

get_beta = function(
  X, # covariates at each site
  x, # observations at each site; can have NAs
  tau, # noise covmats at each sites; dimension may vary ollowing the number of observations
  beta_prior_mean, # prior for the regression coefficients
  beta_prior_precision # prior for the regression coefficients
  )
{
  # Amelioration possible : utiliser la formule de Kronecker pour avoir
  # (chol tau x X) T (chol tau x X) =  (chol tau T x XT)  (chol tau x X) = tau x XTX
  
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
    split(x, row(x)),
    noise_chols))
  # weight each observation of x using noise Choleskys
  chol_noise_x = do.call(c,
    mapply(SIMPLIFY = F,
      function(row_from_x, noise_precision_chol)
    {
      M = matrix(0, length(row_from_x), length(row_from_x))
      M[!is.na(row_from_x), !is.na(row_from_x)] = noise_precision_chol
      row_from_x_ = row_from_x;row_from_x_[is.na(row_from_x_)] = 0
      return(M %*% row_from_x_)
    }, 
    split(x, row(x)),
    noise_chols))
  beta_posterior_precision = beta_prior_precision + crossprod(weighted_X)
  beta_posterior_mean = beta_prior_mean + 
    solve(beta_posterior_precision, crossprod(weighted_X, chol_noise_x))
  return(list("mean" = beta_posterior_mean, "precision" = beta_posterior_precision))
}

# get number of rows of symmetric matrix given number of components
nrow_sqare_matrix = function(mat_coord_length)(-1+ sqrt(1 + 8*mat_coord_length))/2
nrow_sqare_matrix(3)

# get exponential of symmetric matrix
get_symmetric_expmat = function(mat_coordinates)
{
  M = matrix(0,nrow_sqare_matrix(length(mat_coordinates)),nrow_sqare_matrix(length(mat_coordinates)))
  M[lower.tri(M, diag = T)] = mat_coordinates
  M = t(M)
  M[lower.tri(M, diag = T)] = mat_coordinates
  expm::expm(M)
}
  
# get derivatives of exponential of symmetric matrix
get_symmetric_expmat_derivatives = function(mat_coordinates)
{
  M_0 = get_symmetric_expmat(mat_coordinates)
  det_M_0 = determinant(M_0)$mod
  lapply(seq(length(mat_coordinates)), function(i)
  {
    mat_coordinates_ = mat_coordinates; mat_coordinates_[i] = mat_coordinates_[i] + 0.000001
    M = get_symmetric_expmat(mat_coordinates_)
    list("mat_derivative" = 1000000 * (M- M_0), "log_det_derivative" = 1000000 *(determinant(M)$mod - det_M_0))
  }
    )
}

get_symmetric_expmat(c(0,0,0))
get_symmetric_expmat(c(.0001,0,0))
get_symmetric_expmat_derivatives(c(0,0,0))
get_symmetric_expmat_derivatives(c(0))

get_symmetric_expmat_derivatives(c(0,0,0))[[1]][[1]]-(get_symmetric_expmat(c(0.00001,0,0))-get_symmetric_expmat(c(0,0,0)))*100000
get_symmetric_expmat_derivatives(c(0,0,0))[[1]][[2]]-(determinant(get_symmetric_expmat(c(0.00001,0,0)))$mod-determinant(get_symmetric_expmat(c(0,0,0)))$mod)*100000
get_symmetric_expmat_derivatives(c(0,0,0))[[2]][[1]]-(get_symmetric_expmat(c(0,0.00001,0))-get_symmetric_expmat(c(0,0,0)))*100000
get_symmetric_expmat_derivatives(c(0,0,0))[[2]][[2]]-(determinant(get_symmetric_expmat(c(0,0.00001,0)))$mod-determinant(get_symmetric_expmat(c(0,0,0)))$mod)*100000

X_tau = cbind(1, seq(1000)/1000)
beta_tau = matrix(rnorm(6), nrow = 2)
tau_log_coords = X_tau %*% beta_tau

epsilon = t(sapply(split(tau_log_coords, row(tau_log_coords)), function(x)t(chol(solve(get_symmetric_expmat(x))))%*%rnorm(2)))
epsilon[seq(10), 1]=NA
epsilon[seq(11, 20), 2]=NA

epsilon = epsilon + 1 * rnorm(length(epsilon))
# get gradients of beta_tau given a residual epsilon, the current log coordinates of tau, and the regressors X_tau
# uses Jacobian chain rule : 
# nabla beta_tau H = J^T beta tau log tau Nabla log tau H = X_tau Nabla log tau H 
get_beta_tau_gradient = function(epsilon, beta_tau, X_tau)
{
  tau_log_coords = X_tau %*% beta_tau
  # getting NA pattern of the derivative of H wrt log tau coordinates
  NA_pattern =   matrix(t(apply(epsilon, 1, function(epsilon_row)
  {
    res = matrix(1, ncol(epsilon), ncol(epsilon))
    res [is.na(epsilon_row), ] = NA
    res [, is.na(epsilon_row)] = NA
    res[lower.tri(res, diag = T)]
  }
    )), ncol = ncol(epsilon)*(ncol(epsilon)+1)/2)
  derivatives =
    unlist(
      mapply(
        function(tau_log_coord, epsilon_row, NA_pattern_row) 
        {
          lapply(get_symmetric_expmat_derivatives(tau_log_coord[!is.na(NA_pattern_row)]), 
                 function(log_mat_derivative)
                 {
                   return(
                   .5 * sum(epsilon_row[!is.na(epsilon_row)]*(epsilon_row[!is.na(epsilon_row)] %*% log_mat_derivative$mat_derivative)) 
                   - .5 * log_mat_derivative$log_det_derivative
                   )
                 }
                 )
        },
        split(tau_log_coords, row(tau_log_coords)),
        split(epsilon, row(epsilon)), 
        split(NA_pattern, row(NA_pattern))
      )
    )
  derivative_wrt_tau_log_coords = t(NA_pattern)
  derivative_wrt_tau_log_coords[!is.na(derivative_wrt_tau_log_coords)] = derivatives
  derivative_wrt_tau_log_coords[is.na(derivative_wrt_tau_log_coords)] = 0
  t((derivative_wrt_tau_log_coords %*% X_tau))
}

X_tau = cbind(1, seq(1000)/1000)
beta_tau = matrix(rnorm(6), nrow = 2)
tau_log_coords = X_tau %*% beta_tau
epsilon = matrix(sapply(split(tau_log_coords, row(tau_log_coords)), function(x)t(chol(solve(get_symmetric_expmat(x))))%*%rnorm(2)), ncol =2)
beta_tau_ = beta_tau + matrix(c(0, 0, 0, 0, 0, 0.000001), nrow = 2)
1000000*(
get_tau_potential(tau = get_tau(epsilon = epsilon, beta_tau = beta_tau_, X_tau = X_tau), epsilon = epsilon) - 
get_tau_potential(tau = get_tau(epsilon = epsilon, beta_tau = beta_tau, X_tau = X_tau), epsilon = epsilon)
)
get_beta_tau_gradient(epsilon, beta_tau, X_tau)

get_tau = function(epsilon, beta_tau, X_tau)
{
  NA_pattern =   matrix(t(apply(epsilon, 1, function(epsilon_row)
  {
    res = matrix(1, ncol(epsilon), ncol(epsilon))
    res [is.na(epsilon_row), ] = NA
    res [, is.na(epsilon_row)] = NA
    res[lower.tri(res, diag = T)]
  }
  )), ncol = ncol(epsilon)*(ncol(epsilon)+1)/2)
  tau_log_coords = X_tau %*% beta_tau
  mapply(
      function(tau_log_coord, NA_pattern) 
      {
        get_symmetric_expmat(tau_log_coord[!is.na(NA_pattern)])
      },
      split(tau_log_coords, row(tau_log_coords)),
      split(NA_pattern, row(NA_pattern)), 
      SIMPLIFY = F
    )
}

epsilon = matrix(c(1, 1), 1)
beta_tau =matrix(c(1, 1, 0), 1)
X_tau = cbind(1)

.5 * epsilon %*%get_tau(epsilon = epsilon, beta_tau = beta_tau, X_tau = X_tau)[[1]] %*% t(epsilon)- .5* determinant(get_tau(epsilon = epsilon, beta_tau = beta_tau, X_tau = X_tau)[[1]])$mod
get_tau_potential(get_tau(epsilon = epsilon, beta_tau = beta_tau, X_tau = X_tau), epsilon)

get_tau_potential = function(tau, epsilon)
{
  sum(mapply(
    function(epsilon_row, precision_matrix)
    {
      -.5 * determinant(precision_matrix)$mod + .5 * sum(epsilon_row[!is.na(epsilon_row)] * (precision_matrix%*%epsilon_row[!is.na(epsilon_row)]))
    }, 
    split(epsilon, row(epsilon)), 
    tau
  ))
}



# simulating data

locs = seq(1000)
X = cbind(1, seq(1000)/1000)
beta = 3 * matrix(rnorm(4), 2)
X_beta = X%*% beta

X_tau = cbind(1, seq(1000)/1000)
X_tau[,2]= (X_tau[,2]-mean(X_tau[,2]))/sd(X_tau[,2])
beta_tau = matrix(rnorm(6), nrow = 2)+3
tau_log_coords = X_tau %*% beta_tau

epsilon = matrix(t(
  sapply(split(tau_log_coords, row(tau_log_coords)), 
         function(x)t(chol(solve(get_symmetric_expmat(x))))%*%rnorm(2))
  ), ncol = ncol(beta))
observed_x = X_beta + epsilon




state = list()
state$params = list()
state$params$beta = matrix(0, ncol(X), ncol(observed_x))
state$params$beta_tau =  matrix(0, ncol(X_tau), ncol(observed_x)*(ncol(observed_x) + 1)/2)
state$sparse_chol_and_stuff = list()
state$momenta = list()
state$momenta$beta_tau = state$params$beta_tau
state$momenta$beta_tau[] = rnorm(length(state$params$beta_tau))
state$stepsize = list()
state$stepsize$beta_tau = -6
niter = 100
beta_prior_mean = rep(0, ncol(X)*ncol(observed_x))
beta_prior_precision  = diag(1e-6, ncol(X)*ncol(observed_x))
state$sparse_chol_and_stuff$tau = get_tau(epsilon = observed_x, beta_tau = state$params$beta_tau, X_tau = X_tau)


for(iter in seq(niter)){
   print(iter) 
  # updating noise
  proposed_beta_tau = state$params$beta_tau
  print(state$params$beta_tau)
  state$momenta$beta_tau = state$momenta$beta_tau * 
    sqrt(.9) + rnorm(length(beta_tau)) * sqrt(.1)
  proposed_beta_tau_momentum = state$momenta$beta_tau
  
  current_tau_potential = get_tau_potential(tau = get_tau(epsilon = observed_x, beta_tau = proposed_beta_tau, X_tau = X_tau), epsilon = observed_x - X%*%beta)
  current_tau_kinetic = sum(state$momenta$beta_tau^2)/2
  t1 = Sys.time()
  #half step for momentum
  proposed_beta_tau_momentum = proposed_beta_tau_momentum - .5 * exp(state$stepsize$beta_tau) * get_beta_tau_gradient(epsilon = observed_x - X%*%beta, beta_tau = proposed_beta_tau, X_tau = X_tau)
  #full step to update beta tau
  proposed_beta_tau = proposed_beta_tau + exp(state$stepsize$beta_tau) * proposed_beta_tau_momentum
  #half step for momentum
  proposed_beta_tau_momentum = proposed_beta_tau_momentum - .5 * exp(state$stepsize$beta_tau) * get_beta_tau_gradient(epsilon = observed_x - X%*%beta, beta_tau = proposed_beta_tau, X_tau = X_tau)
  print(Sys.time()-t1)
  
  proposed_tau_potential = get_tau_potential(tau = get_tau(epsilon = observed_x, beta_tau = proposed_beta_tau, X_tau = X_tau), epsilon = observed_x - X%*%beta)
  proposed_tau_kinetic = sum(proposed_beta_tau_momentum^2)/2
  
  
  if(- proposed_tau_potential - proposed_tau_kinetic + current_tau_potential + current_tau_kinetic > log(runif(1)))
  {
    print("tatatooooo !")
    state$momenta$beta_tau = proposed_beta_tau_momentum
    state$params$beta_tau = proposed_beta_tau
  }
    
  
  state$sparse_chol_and_stuff$tau = get_tau(epsilon = observed_x, beta_tau = state$params$beta_tau, X_tau = X_tau)
  
  betas_distribution = get_beta (
    X = X, # covariates at each site
    x = observed_x, # observations at each site; can have NAs
    tau = state$sparse_chol_and_stuff$tau, # noise covmats at each sites; dimension may vary ollowing the number of observations
    beta_prior_mean = beta_prior_mean, # prior for the regression coefficients
    beta_prior_precision  = beta_prior_precision # prior for the regression coefficients
  )
  state$params$beta[] = betas_distribution$mean + solve(chol(betas_distribution$precision), rnorm(length(betas_distribution$mean)))
}





M = matrix(0, 3, 3)
tatato = rnorm(6)
M[lower.tri(M, T)] = tatato; M = t(M); M[lower.tri(M, T)] = tatato
eig = eigen(M)
eig$vectors%*%diag(eig$values)%*%t(eig$vectors)-M

tcrossprod(eig$vectors[,1])*eig$values[1]+
tcrossprod(eig$vectors[,2])*eig$values[2]+
tcrossprod(eig$vectors[,3])*eig$values[3]
jac = do.call(rbind, list(
  tcrossprod(eig$vectors[,1])[lower.tri(tcrossprod(eig$vectors[,1]), T)], 
  tcrossprod(eig$vectors[,2])[lower.tri(tcrossprod(eig$vectors[,1]), T)], 
  tcrossprod(eig$vectors[,3])[lower.tri(tcrossprod(eig$vectors[,1]), T)],
  (eig$values[1] * (tcrossprod(eig$vectors[,1], c(1, 0, 0)) + tcrossprod(c(1, 0, 0), eig$vectors[,1]))[lower.tri(tcrossprod(eig$vectors[,1]), T)]),
  (eig$values[1] * (tcrossprod(eig$vectors[,1], c(0, 1, 0)) + tcrossprod(c(0, 1, 0), eig$vectors[,1]))[lower.tri(tcrossprod(eig$vectors[,1]), T)])
  )
)
M_ = M
i = 1
j = 1
M_[i,j]= M_[i,j]+.00001
eig_ = eigen(M_)
((eig_$vectors-eig$vectors)*100000)
(eig_$values-eig$values)*100000
jac
pracma::pinv(jac)


eig__ = eig
eig__$values[1]= eig__$values[1]+0.001
(eig__$vectors%*%diag(eig__$values)%*%t(eig__$vectors)-M)*1000

eig__ = eig
eig__$vectors[1, 1]= eig__$vectors[1, 1]+0.001
(eig__$vectors%*%diag(eig__$values)%*%t(eig__$vectors)-M)*1000
jac

eig__ = eig
eig__$vectors[2, 1]= eig__$vectors[2, 1]+0.001
(eig__$vectors%*%diag(eig__$values)%*%t(eig__$vectors)-M)*1000
jac




M = matrix(0, 3, 3)
tatato = rnorm(6)
M[lower.tri(M, T)] = tatato; M = t(M); M[lower.tri(M, T)] = tatato
eig = eigen(M)
eig$vectors%*%diag(eig$values)%*%t(eig$vectors)-M

i = 1
j = 2 
k = 3
dMdMjk = matrix(0, nrow(M), nrow(M));dMdMjk[k, j]=1;dMdMjk[j,k]=1
D_ii_derivative = t(eig$vectors[,i]) %*% dMdMjk %*% eig$vectors[,i]

M_ = M; M_[j, k] =M_[j, k]+0.001; M_[k, j] =M_[j, k]
eig_ = eigen(M_)
print(D_ii_derivative)
print((eig_$values[i]-eig$values[i])*1000)

pseudoinv_diag = 1/eig$values; pseudoinv_diag[i]= 0.000000
pseudoinv = eig$vectors%*%diag(pseudoinv_diag)%*%t(eig$vectors)

pseudoinv %*% (dMdMjk - c(D_ii_derivative) * diag(1, nrow(M), nrow(M))) %*% eig$vectors[,i]
(eig_$vectors[,i]-eig$vectors[,i])*1000

-(dMdMjk - c(D_ii_derivative) * diag(1, nrow(M), nrow(M))) %*% eig$vectors[,i]
(M - eig$values[i] * diag(1, nrow(M), nrow(M)))%*%(eig_$vectors[,i]-eig$vectors[,i])*1000

- solve(M - eig$values[i] * diag(1, nrow(M), nrow(M)) + diag(0.00000001, nrow(M), nrow(M))) %*% (dMdMjk - c(D_ii_derivative) * diag(1, nrow(M), nrow(M))) %*% eig$vectors[,i]
(eig_$vectors[,i]-eig$vectors[,i])*1000

eigen(M - eig$values[i] * diag(1, nrow(M), nrow(M)))
solve(M - eig$values[i] * diag(1, nrow(M), nrow(M)) + diag(0.00000001, nrow(M), nrow(M)))
1/eig$values

(M - eig$values[i] * diag(1, nrow(M), nrow(M))) %*% pseudoinv





