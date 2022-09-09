set.seed(1)

#space_locs = cbind(runif(100), runif(100))
space_locs = cbind(rep(runif(10), 10), rep(runif(10), 10))
time_locs = cbind(seq(10))
p = 3
var_tag = rep(1+floor(p*runif(10)), 10)
a2 = runif(p)
b2 = runif(p)
nu = .5 + c(0, 1, 2.5)
nu_ = outer(nu, nu, "+"); nu_ = nu_/2; nu_ = nu_[lower.tri(nu_, diag = T)]
alpha = .0001 * runif(1)
beta = runif(1)
rho = crossprod(matrix(rnorm(p^2), p))
gamma_0_a = runif(1) 
gamma_0_b = runif(1) 
gamma_0_c = 50 * runif(1) 
R_lambda = runif(1) 
R_b = runif(1) 
R_r = 50 * runif(1) 
A = runif(p)
u = seq(0, 9)

eta_and_derivatives = function(u, # unique vector of time distance
                               A, # vector of size p with elements between 0 and 1, controlling time correlation
                               gamma_0_c, # number greater than 0, controlling time correlation
                               gamma_0_a, gamma_0_b, # numbers between 0 and 1, controlling time correlation
                               R_c, R_lambda, R_b # numbers between 0 and 1, controlling time correlation
)
  {
  gamma_0 = ((1 + (gamma_0_c*u)^(2*gamma_0_a))^gamma_0_b)
  gamma_0_deriv = eval(deriv(~((1 + (gamma_0_c*u)^(2*gamma_0_a))^gamma_0_b), name = c("gamma_0_c", "gamma_0_a", "gamma_0_b"), hessian = T))
  time_effect = ((1 + (R_r * u)^(2*R_lambda))^(-R_b))
  time_effect_deriv = eval(deriv(~((1 + (R_r * u)^(2*R_lambda))^(-R_b)), name = c("R_r", "R_lambda", "R_b"), hessian = T))
  A_mat = outer(A, A); vech_A_mat = A_mat[lower.tri(A_mat, diag = T)]
  eta = 
  outer(gamma_0, rep(1, length(vech_A_mat)))+
  outer(rep(1, length(u)), vech_A_mat)+
  outer(time_effect, vech_A_mat)
  # gradient of eta:
  # row-constant wrt gamma 0 's parameters
  # proportional to vech_A_mat wrt time effect 's parameters
  # hessian of eta : 
  # null when crossing gamma 0 and time effects 's parameters
  # null when crossing gamma 0 and A
  # When derivative wrt A: use sparsity !
  return(list("eta" = eta, "vech_A_mat" = vech_A_mat, "time_effect" = time_effect, "time_effect_deriv" = time_effect_deriv, "gamma_0" = gamma_0, "gamma_0_deriv" = gamma_0_deriv))
}

compute_eta = function(u, # unique vector of time distance
                               A, # vector of size p with elements between 0 and 1, controlling time correlation
                               gamma_0_c, # number greater than 0, controlling time correlation
                               gamma_0_a, gamma_0_b, # numbers between 0 and 1, controlling time correlation
                               R_c, R_lambda, R_b # numbers between 0 and 1, controlling time correlation
)
  {
  gamma_0 = ((1 + (gamma_0_c*u)^(2*gamma_0_a))^gamma_0_b)
  time_effect = (1- # R0
                   (1 + (R_r * u)^(2*R_lambda))^(-R_b) # -R(u)
                 )
  A_mat = outer(A, A); vech_A_mat = A_mat[lower.tri(A_mat, diag = T)]
  eta = 
  outer(gamma_0, rep(1, length(vech_A_mat)))+
  outer(time_effect, vech_A_mat)
  dimnames(eta)[[1]] = paste("u=", u, sep = "_")
  dimnames(eta)[[2]] = paste("v1=",unlist(sapply(seq(length(A)), function(x) seq(x, length(A)))), "_v2=",  rep(seq(length(A)), seq(length(A), 1)), sep = "_")
  return(eta)
}

compute_eta_and_variation = function(u, # unique vector of time distance
                               A, # vector of size p with elements between 0 and 1, controlling time correlation
                               gamma_0_c, # number greater than 0, controlling time correlation
                               gamma_0_a, gamma_0_b, # numbers between 0 and 1, controlling time correlation
                               R_c, R_lambda, R_b # numbers between 0 and 1, controlling time correlation
)
  {
  eta = compute_eta(u = u, A = A, gamma_0_c = gamma_0_c, gamma_0_a = gamma_0_a, gamma_0_b = gamma_0_b, R_c = R_c, R_lambda = R_lambda, R_b = R_b)
  eta_variation = array(0, dim = c(dim(eta), length(A)+6))
  eta_variation[,,1] = compute_eta(u = u, A = A, gamma_0_c = gamma_0_c + .0001, gamma_0_a = gamma_0_a, gamma_0_b = gamma_0_b, R_c = R_c, R_lambda = R_lambda, R_b = R_b)
  eta_variation[,,2] = compute_eta(u = u, A = A, gamma_0_c = gamma_0_c, gamma_0_a = gamma_0_a + .0001, gamma_0_b = gamma_0_b, R_c = R_c, R_lambda = R_lambda, R_b = R_b)
  eta_variation[,,3] = compute_eta(u = u, A = A, gamma_0_c = gamma_0_c, gamma_0_a = gamma_0_a, gamma_0_b = gamma_0_b + .0001, R_c = R_c, R_lambda = R_lambda, R_b = R_b)
  eta_variation[,,4] = compute_eta(u = u, A = A, gamma_0_c = gamma_0_c, gamma_0_a = gamma_0_a, gamma_0_b = gamma_0_b, R_c = R_c + .0001, R_lambda = R_lambda, R_b = R_b)
  eta_variation[,,5] = compute_eta(u = u, A = A, gamma_0_c = gamma_0_c, gamma_0_a = gamma_0_a, gamma_0_b = gamma_0_b, R_c = R_c, R_lambda = R_lambda + .0001, R_b = R_b)
  eta_variation[,,6] = compute_eta(u = u, A = A, gamma_0_c = gamma_0_c, gamma_0_a = gamma_0_a, gamma_0_b = gamma_0_b, R_c = R_c, R_lambda = R_lambda, R_b = R_b + .0001)
  for(i in seq(length(A)))eta_variation[,,6] =
    {
      A_ = A
      A_[i] = A[i] + .0001
      eta_variation[,,i + 6] =compute_eta(u = u, A = A_, gamma_0_c = gamma_0_c, gamma_0_a = gamma_0_a, gamma_0_b = gamma_0_b, R_c = R_c, R_lambda = R_lambda, R_b = R_b)
    }
  dimnames(eta_variation)[[1]] = paste("u=", u, sep = "_")
  dimnames(eta_variation)[[2]] = paste("v1=",unlist(sapply(seq(length(A)), function(x) seq(x, length(A)))), "_v2=",  rep(seq(length(A)), seq(length(A), 1)), sep = "_")
  dimnames(eta_variation)[[3]] = c("gamma_0_c", "gamma_0_a", "gamma_0_b", "R_c", "R_lambda", "R_b", paste("A", seq(length(A)), sep = "_"))
  return(list("eta" = eta, "eta_variation" = eta_variation))
}


eta_and_variation = compute_eta_and_variation(u = u, A = A, 
                                              gamma_0_c = gamma_0_c, gamma_0_a = gamma_0_a, gamma_0_b = gamma_0_b, R_c = R_c, R_lambda = R_lambda, R_b = R_b )
eta = eta_and_variation[[1]]
compute_R = function(eta, a2, b2, alpha, beta)
{
  # expanding range parameters a2 and b2 
  a2_ = outer(a2, a2, "+"); a2_ = a2_/2; a2_ = a2_[lower.tri(a2_, diag = T)]
  b2_ = outer(b2, b2, "+"); b2_ = b2_/2; b2_ = b2_[lower.tri(b2_, diag = T)]
  R = sqrt((alpha * eta + matrix(1, nrow(eta), 1) %*% a2_)/(beta * eta + matrix(1, nrow(eta), 1) %*% b2_ ))
}


compute_R_and_variation = function(eta_and_variation, a2, b2, alpha, beta){
  R = compute_R(eta_and_variation[[1]], a2, b2, alpha, beta)
  R_variation = array(0, dim = c(dim(eta_and_variation[[1]]),
  dim(eta_and_variation$eta_variation)[3]+2*(length(a2)+1) + length(a2)))
  R_variation[,,seq(length(a2)+6)] = apply(
    eta_and_variation$eta_variation, 3, function(x)compute_R(x, a2, b2, alpha, beta)
  )
  R_variation[,,length(a2)+7] = compute_R(eta_and_variation[[1]], a2, b2, alpha + 0.0001, beta)
  R_variation[,,length(a2)+8] = compute_R(eta_and_variation[[1]], a2, b2, alpha, beta + 0.0001)
  for(i in seq(length(a2))){
    a2_ = a2; a2_[i] = a2_[i] + .0001
    b2_ = b2; b2_[i] = b2_[i] + .0001
    R_variation[,,  length(a2)+8+i] = compute_R(eta_and_variation[[1]], a2_, b2, alpha, beta)
    R_variation[,,2*length(a2)+8+i] = compute_R(eta_and_variation[[1]], a2, b2_, alpha, beta)
  }
  for(i in seq(length(a2))){
    R_variation[,,3*length(a2)+8+i] = R
  }
  dimnames(R_variation)[[1]] = dimnames(eta_and_variation$eta_variation)[[1]]
  dimnames(R_variation)[[2]] = dimnames(eta_and_variation$eta_variation)[[2]]
  dimnames(R_variation)[[3]] = c(dimnames(eta_and_variation$eta_variation)[[3]], "alpha", "beta", paste("a2", seq(length(a2)), sep = "_"), paste("b2", seq(length(a2)), sep = "_"), paste("nu", seq(length(a2)), sep = "_"))
  return(list("R" = R, "R_variation" = R_variation))
}



R_and_variation = compute_R_and_variation(eta_and_variation, a2, b2, alpha, beta)


Matern = function(h, r, nu_)(2^(1-nu_)) * (r * h)^nu_ * besselK(r * h, nu_)

compute_sig = function(eta, a2, b2, alpha, beta, nu)
{
  # expanding range parameters a2 and b2 
  a2_ = outer(a2, a2, "+"); a2_ = a2_/2; a2_ = a2_[lower.tri(a2_, diag = T)]
  b2_ = outer(b2, b2, "+"); b2_ = b2_/2; b2_ = b2_[lower.tri(b2_, diag = T)]
  nu_ = outer(nu, nu, "+"); nu_ = nu_/2; nu_ = nu_[lower.tri(nu_, diag = T)]
  sig_normalization =  
    (1/gamma(nu)) * # normalizing diag of gamma(nu) for both matern and sig. When nu_i \approx nu_j, close to 1. 
    sqrt((alpha + a2)^nu) * #diag(eta) = 1
    sqrt((beta + b2)^(ncol(space_locs)/2)) #diag(eta) = 1
  sig_normalization = outer(sig_normalization, sig_normalization)
  sig_normalization = sig_normalization[lower.tri(sig_normalization, diag = T)]
  sig = 1/
    (
      (alpha * eta + outer(rep(1, nrow(eta)),a2_) )^outer(rep(1, nrow(eta)),nu_) *
        (beta * eta + outer(rep(1, nrow(eta)),b2_))^(ncol(space_locs)/2)
    )*
    outer(rep(1, nrow(eta)),gamma(nu_)) * 
    # * exp(-nu) removed because of negative definiteness
  outer(rep(1, nrow(eta)), sig_normalization)
  sig
}

compute_sig_and_variation = function(eta_and_variation, a2, b2, alpha, beta, nu)
{
  sig = compute_sig(eta = eta_and_variation[[1]], a2 = a2, b2 = b2, alpha = alpha, beta = beta, nu = nu)
  sig_variation = array(0, dim = c(dim(eta_and_variation[[1]]),
                                 dim(eta_and_variation$eta_variation)[3]+3*(length(a2))+2))
  
  sig_variation[,,seq(length(a2)+6)] = apply(
    eta_and_variation$eta_variation, 3, function(x)compute_sig(x, a2, b2, alpha, beta, nu)
  )
  sig_variation[,,length(a2)+7] = compute_sig(eta_and_variation[[1]], a2, b2, alpha + 0.0001, beta, nu)
  sig_variation[,,length(a2)+8] = compute_sig(eta_and_variation[[1]], a2, b2, alpha, beta + 0.0001, nu)
  for(i in seq(length(a2))){
    a2_ = a2; a2_[i] = a2_[i] + .0001
    b2_ = b2; b2_[i] = b2_[i] + .0001
    nu_ = nu; nu_[i] = nu_[i] + .0001
    sig_variation[,,  length(a2)+8+i] = compute_sig(eta_and_variation[[1]], a2_, b2, alpha, beta, nu)
    sig_variation[,,2*length(a2)+8+i] = compute_sig(eta_and_variation[[1]], a2, b2_, alpha, beta, nu)
    sig_variation[,,3*length(a2)+8+i] = compute_sig(eta_and_variation[[1]], a2, b2, alpha, beta, nu_)
  }
  dimnames(sig_variation)[[1]] = dimnames(eta_and_variation$eta_variation)[[1]]
  dimnames(sig_variation)[[2]] = dimnames(eta_and_variation$eta_variation)[[2]]
  dimnames(sig_variation)[[3]] = c(dimnames(eta_and_variation$eta_variation)[[3]], "alpha", "beta", paste("a2", seq(length(a2)), sep = "_"), paste("b2", seq(length(a2)), sep = "_"), paste("nu", seq(length(a2)), sep = "_"))
  return(list("sig" = sig, "sig_variation" = sig_variation))
}  

sig_and_variation = compute_sig_and_variation(eta_and_variation, a2, b2, alpha, beta, nu)

# given integer n, gives 2-column matrix of lower triangular indices 
# in the square matrix n*n WITHOUT diagonal as in dist(,diag = F)
lower_tri_idx = function(n, diag = F)
{
  
  if(diag == F)
  {
    return(
      cbind(
       rev(abs(sequence(seq.int(n - 1)) - n) + 1),
       rep.int(seq.int(n - 1), rev(seq.int(n - 1)))
      )
    )
  }
  if(diag == T)
  {
    return(
      cbind(
       rev(abs(sequence(seq.int(n)) - n) + 1),
       rep.int(seq.int(n), rev(seq.int(n)))
      )
    )
  }
}
lower_tri_idx(10)

# given three integers i, j, n, gives the position of coefficient (i, j)
# in the n*n square matrix in the lower triangular coefficients
# WITH diagonal as in dist(,diag = T)
position_in_lower_tri = function(i, j, n)
{  
  a = pmin(i, j); b = pmax(i, j)
  (a-1)*(n-a+1) + (a-1)*a/2 + b-a+1
}
i = 1;j =3;n=3;position_in_lower_tri(i, j, n)
i = 2;j =1;n=3;position_in_lower_tri(i, j, n)
i = 2;j =3;n=3;position_in_lower_tri(i, j, n)
i = 3;j =3;n=3;position_in_lower_tri(i, j, n)


# given a vector of integers i_vec and an integer n_var, 
# n var being greater than the entries of i_vec, 
# computes the lower-triangular without diagonal of the crossing matrix of i_vec, i_vec
# and finds the position of those pairs in the lower triangular part of the 
# matrix of size n_var * n_var
position_in_lower_tri_cross_vec = function(i_vec, n_var, diag = F)
{  
  idx = lower_tri_idx(length(i_vec), diag = diag)
  position_in_lower_tri(i_vec[idx[,1]], i_vec[idx[,2]], n_var)
}

n_var = 3
i_vec = 1 + floor(3*runif(10))
position_in_lower_tri_cross_vec(i_vec, 3)
length(position_in_lower_tri_cross_vec(i_vec, 3))


children = cbind(runif(3), runif(3))
children_var_tag = 1 + floor(n_var*runif(3))
parents_current_time = cbind(runif(80), runif(80))
parents_current_time_var_tag = 1 + floor(n_var*runif(80))
parents_previous_times = cbind(runif(80), runif(80))
parents_previous_times_var_tag = 1 + floor(n_var*runif(80))


rho = crossprod(matrix(rnorm(n_var^2), n_var))
rho = diag(1/sqrt(diag(rho))) %*% rho %*% diag(1/sqrt(diag(rho)))
rho = rho[lower.tri(rho, diag = T)] 


R = R_and_variation[[1]]
sig = sig_and_variation[[1]]

GM = function(
    children, children_var_tag,  # list of children conditionned by parents at the same time and in previous times
    parents_current_time, parents_current_time_var_tag, # list of parents at the same time
    parents_previous_times, parents_previous_times_var_tag, # list of parents at previous times
    R, sig, rho, nu, 
    n_var
)
{
  #getting cross-variable smoothness from margnial smoothness parameters
  nu_ = outer(nu, nu, "+"); nu_ = nu_/2; nu_ = nu_[lower.tri(nu_, diag = T)]
  
  # size of spatial sets of interest
  n1 = nrow(children)
  n2 = nrow(parents_current_time)
  n3 = nrow(parents_previous_times)
  
  # allocating covmat
  covmat = matrix(0, n1+n2+n3*(nrow(R)-1), n1+n2+n3*(nrow(R)-1))
  
  
  # filling covmat for children when there is more than one child 
  if(nrow(children)>1){
    # index of variable pairs in the list of variable pairs 
    var_idx = position_in_lower_tri_cross_vec(children_var_tag, n_var)
    h_children = c(dist(children))
    # multivariate nonseparable GM
    covmat[lower_tri_idx(n1)] = 
    Matern(h_children, r = R[1,var_idx], nu_ = nu_[var_idx]) * 
      sig[1,var_idx] * rho[var_idx]
    # filling transpose 
    covmat[lower_tri_idx(n1)[,c(2, 1), drop = F]] = covmat[lower_tri_idx(n1)]
  }
  covmat[cbind(seq(n1), seq(n1))] = 1
  # checking positive definiteness
  # chol(covmat[seq(n1), seq(n1)])

  
  # parents at same time
  var_idx = position_in_lower_tri_cross_vec(parents_current_time_var_tag, n_var)
  h = c(dist(parents_current_time))
  covmat[n1+lower_tri_idx(n2)] = 
    Matern(h, r = R[1,var_idx], nu_ = nu_[var_idx]) * 
    sig[1,var_idx] * rho[var_idx]
  covmat[n1+lower_tri_idx(n2)[,c(2, 1)]] = covmat[n1+lower_tri_idx(n2)]
  covmat[n1+cbind(seq(n2), seq(n2))] = 1
  #chol(covmat[seq(n1+1, n2), seq(n1+1, n2)])
  
  # parents at previous times
  # var combinations and distances between unique parent pairs
  var_idx = position_in_lower_tri_cross_vec(parents_previous_times_var_tag, n_var, diag = T)
  h = as.matrix(dist(parents_previous_times, diag = T));h = h[lower.tri(h, diag = T)]; 
  h[h==0] = min(h[h!=0])*.00001
  idx = n1+n2+lower_tri_idx(n3, diag = T);idx[,1] = idx[,1]-n3
  for(i in seq(length(u)-1)){
    # index of position in the covariance matrix
    h[h==0] = min(h[h!=0])*.00001
    idx[,1] = idx[,1] + n3
    covmat[idx] = 
      Matern(h, r = R[i,var_idx], nu_ = nu_[var_idx]) * 
      sig[i,var_idx] * rho[var_idx]
    idx_ = n1+n2+lower_tri_idx(n3, diag = T)[,c(2, 1)]; idx_[,1] = idx_[,1] + n3*(i-1)
    covmat[idx_] = covmat[idx]
    covmat[seq(n1+n2+1,n1+n2+n3), seq(n1+n2+1+n3*(i-1), n1+n2+n3*(i))] = covmat[seq(n1+n2+1+n3*(i-1), n1+n2+n3*(i)) , seq(n1+n2+1,n1+n2+n3)] 
    if(i!=(length(u)-1))for(j in seq(length(u)-1-i))
    {
      covmat[j*n3 + seq(n1+n2+1+n3*(i-1), n1+n2+n3*(i)) , j*n3 +  seq(n1+n2+1,n1+n2+n3)] =   covmat[seq(n1+n2+1+n3*(i-1), n1+n2+n3*(i)) , seq(n1+n2+1,n1+n2+n3)]
      covmat[j*n3 +  seq(n1+n2+1,n1+n2+n3) , j*n3 + seq(n1+n2+1+n3*(i-1), n1+n2+n3*(i))] =   covmat[seq(n1+n2+1+n3*(i-1), n1+n2+n3*(i)) , seq(n1+n2+1,n1+n2+n3)]
    }
  }
  #chol(covmat[-seq(n1+n2), -seq(n1+n2)])
    
  # children * parents current time 
  h = fields::rdist(children, parents_current_time)
  h[h==0] = min(h[h!=0])*.00001
  var_idx = position_in_lower_tri(outer(children_var_tag, rep(1, n2)), outer(rep(1, n1), c(parents_current_time_var_tag)), n_var)
  covmat[seq(n1), seq(n1+1, n1+n2)] = 
    Matern(h, r = matrix(R[1, var_idx], n1), nu_ = matrix(nu_[var_idx], n1)) * 
    sig[1, var_idx]  * rho[var_idx]
  covmat[seq(n1+1, n1+n2), seq(n1)] = t(covmat[seq(n1), seq(n1+1, n1+n2)])

  # children * parents previous time 
  h = fields::rdist(children, parents_previous_times)
  h[h==0] = min(h[h!=0])*.00001
  var_idx = position_in_lower_tri(outer(children_var_tag, rep(1, n3)), outer(rep(1, n1), c(parents_previous_times_var_tag)), n_var)
  for(i in seq(2, length(u)))
  {
    covmat[seq(n1), seq(n1+n2+(i-2)*n3+1, n1+n2+(i-1)*n3)] = 
      Matern(h, r = matrix(R[i, var_idx], n1), nu_ = matrix(nu_[var_idx], n1)) * 
      sig[i, var_idx]  * rho[var_idx]
    covmat[seq(n1+n2+(i-2)*n3+1, n1+n2+(i-1)*n3), seq(n1)] = covmat[seq(n1), seq(n1+n2+(i-2)*n3+1, n1+n2+(i-1)*n3)]
  }

  # parents current time * parents previous time 
  h = fields::rdist(parents_current_time, parents_previous_times)
  h[h==0] = min(h[h!=0])*.00001
  var_idx = position_in_lower_tri(outer(parents_current_time_var_tag, rep(1, n3)), outer(rep(1, n2), c(parents_previous_times_var_tag)), n_var)
  for(i in seq(2, length(u)))
  {
    covmat[seq(n1+1,n1+n2), seq(n1+n2+(i-2)*n3+1, n1+n2+(i-1)*n3)] = 
      Matern(h, r = matrix(R[i, var_idx], n2), nu_ = matrix(nu_[var_idx], n2)) * 
      sig[i, var_idx]  * rho[var_idx]
    covmat[seq(n1+n2+(i-2)*n3+1, n1+n2+(i-1)*n3) , seq(n1+1,n1+n2)] = t(covmat[seq(n1+1,n1+n2), seq(n1+n2+(i-2)*n3+1, n1+n2+(i-1)*n3)])
  }
  covmat
  #chol(covmat)
  #image(covmat)
}

t1 = Sys.time()
for(i in seq(100))
{
  tatato =
GM(children = children, children_var_tag = children_var_tag, 
   parents_current_time = parents_current_time, parents_current_time_var_tag = parents_current_time_var_tag, 
   parents_previous_times = parents_previous_times, parents_previous_times_var_tag = parents_previous_times_var_tag, 
   R = R, sig = sig, rho = rho, nu = nu, n_var = n_var)
}
Sys.time()-t1

t1 = Sys.time()
for(i in seq(100))
{
  
tatato = cbind(rnorm(883), rnorm(883))
tatato = GpGp::matern_isotropic(c(1, 1, 1, 0), tatato)
}
Sys.time()-t1



GM_ = function(
    children, children_var_tag,  # list of children conditionned by parents at the same time and in previous times
    parents_current_time, parents_current_time_var_tag, # list of parents at the same time
    parents_previous_times, parents_previous_times_var_tag, # list of parents at previous times
    R_and_variation, sig_and_variation,
    rho, nu, 
    n_var
)
{
  
  #getting cross-variable smoothness from margnial smoothness parameters
  nu_ = outer(nu, nu, "+"); nu_ = nu_/2; nu_ = nu_[lower.tri(nu_, diag = T)]
  
  # size of spatial sets of interest
  n1 = nrow(children)
  n2 = nrow(parents_current_time)
  n3 = nrow(parents_previous_times)
  
  res = list()
  res$value = matrix(0, n1, n1+n1+n3*(length(u)-1))
  res$derivative = array(0, dim =c(n1, n1+n1+n3*(length(u)-1), 4*n_var+n_var*(n_var-1)/2))
  dimnames(res$derivative) = list(NULL, NULL, c(
    paste("nu", seq(n_var), sep = "_"), 
    paste("a2", seq(n_var), sep = "_"), 
    paste("b2", seq(n_var), sep = "_"), 
    paste("A", seq(n_var), sep = "_"), 
    paste("rho", outer(seq(n_var), seq(n_var), paste, sep = "")[lower.tri(outer(seq(n_var), seq(n_var), paste, sep = ""))], sep = "_")
    ))
  
  # allocating covmat
  
  covmat_c = matrix(0, n1, n1)
  covmat_p = matrix(0, n2+n3*(length(u)-1), n2+n3*(length(u)-1))
  covmat_cp = matrix(0, n1, n2+n3*(length(u)-1))
  

  covmat_cp_derivative = array(0, dim =c(n1, n2+n3*(length(u)-1), 4*n_var+n_var*(n_var-1)/2))
  dimnames(covmat_cp_derivative) = list(NULL, NULL, c(
    paste("nu", seq(n_var), sep = "_"), 
    paste("a2", seq(n_var), sep = "_"), 
    paste("b2", seq(n_var), sep = "_"), 
    paste("A", seq(n_var), sep = "_"), 
    paste("rho", outer(seq(n_var), seq(n_var), paste, sep = "")[lower.tri(outer(seq(n_var), seq(n_var), paste, sep = ""))], sep = "_")
    ))
  covmat_p_derivative = array(0, dim =c(n2+n3*(length(u)-1), n2+n3*(length(u)-1), 4*n_var+n_var*(n_var-1)/2))
  dimnames(covmat_p_derivative) = list(NULL, NULL, c(
    paste("nu", seq(n_var), sep = "_"), 
    paste("a2", seq(n_var), sep = "_"), 
    paste("b2", seq(n_var), sep = "_"), 
    paste("A", seq(n_var), sep = "_"), 
    paste("rho", outer(seq(n_var), seq(n_var), paste, sep = "")[lower.tri(outer(seq(n_var), seq(n_var), paste, sep = ""))], sep = "_")
    ))
  
  
  
  # filling covmat for children when there is more than one child 
  if(nrow(children)>1){
    # index of variable pairs in the list of variable pairs 
    var_idx = position_in_lower_tri_cross_vec(children_var_tag, n_var)
    h = c(dist(children))
    h[h==0] = min(h[h!=0])*.00001
    # multivariate nonseparable GM
    covmat_c[lower_tri_idx(n1)] = 
      Matern(h, r = R[1,var_idx], nu_ = nu_[var_idx]) * 
      sig[1,var_idx] * rho[var_idx]
    # filling transpose 
    covmat_c[lower_tri_idx(n1)[,c(2, 1), drop = F]] = covmat_c[lower_tri_idx(n1)]
  }
  covmat_c[cbind(seq(n1), seq(n1))] = 1
  # checking positive definiteness
  # chol(covmat_c)
  
  
  # children * parents current time 
  h = fields::rdist(children, parents_current_time)
  h[h==0] = min(h[h!=0])*.00001
  var_idx = position_in_lower_tri(outer(children_var_tag, rep(1, n2)), outer(rep(1, n1), c(parents_current_time_var_tag)), n_var)
  covmat_cp[seq(n1), seq(n2)] = 
    Matern(h, r = matrix(R[1, var_idx], n1), nu_ = matrix(nu_[var_idx], n1)) * 
    sig[1, var_idx]  * rho[var_idx]
  
  # children * parents previous time 
  h = fields::rdist(children, parents_previous_times)
  h[h==0] = min(h[h!=0])*.00001
  var_idx = position_in_lower_tri(outer(children_var_tag, rep(1, n3)), outer(rep(1, n1), c(parents_previous_times_var_tag)), n_var)
  for(i in seq(2, length(u)))
  {
    covmat_cp[seq(n1), seq(n2+(i-2)*n3+1,n2+(i-1)*n3)] = 
      Matern(h, r = matrix(R[i, var_idx], n1), nu_ = matrix(nu_[var_idx], n1)) * 
      sig[i, var_idx]  * rho[var_idx]
  }
  #image(covmat_cp)
  
  
  # parents at same time
  var_idx = position_in_lower_tri_cross_vec(parents_current_time_var_tag, n_var)
  h = c(dist(parents_current_time))
  h[h==0] = min(h[h!=0])*.00001
  mat_idx = lower_tri_idx(n2)
  covmat_p[mat_idx] = 
    Matern(h, r = R[1,var_idx], nu_ = nu_[var_idx]) * 
    sig[1,var_idx] * rho[var_idx]
  covmat_p[mat_idx[,c(2, 1)]] = covmat_p[mat_idx]
  covmat_p[cbind(seq(n2), seq(n2))] = 1
  #chol(covmat_p[seq(n2), seq(n2)])
  
  
  
  
  # parents at previous times
  # var combinations and distances between unique parent pairs
  var_idx = position_in_lower_tri_cross_vec(parents_previous_times_var_tag, n_var, diag = T)
  h = as.matrix(dist(parents_previous_times, diag = T));h = h[lower.tri(h, diag = T)]; 
  h[h==0] = min(h[h!=0])*.00001
  idx = n2+lower_tri_idx(n3, diag = T);idx[,1] = idx[,1]-n3
  for(i in seq(length(u)-1)){
    # index of position in the covariance matrix
    h[h==0] = min(h[h!=0])*.00001
    idx[,1] = idx[,1] + n3
    covmat_p[idx] = 
      Matern(h, r = R[i,var_idx], nu_ = nu_[var_idx]) * 
      sig[i,var_idx] * rho[var_idx]
    idx_ = n2+lower_tri_idx(n3, diag = T)[,c(2, 1)]; idx_[,1] = idx_[,1] + n3*(i-1)
    covmat_p[idx_] = covmat_p[idx]
    covmat_p[seq(n2+1,n2+n3), seq(n2+1+n3*(i-1), n2+n3*(i))] = covmat_p[seq(n2+1+n3*(i-1), n2+n3*(i)) , seq(n2+1,n2+n3)] 
    if(i!=(length(u)-1))for(j in seq(length(u)-1-i))
    {
      covmat_p[j*n3 + seq(n2+1+n3*(i-1), n2+n3*(i)) , j*n3 +  seq(n2+1,n2+n3)] =   covmat_p[seq(n2+1+n3*(i-1), n2+n3*(i)) , seq(n2+1,n2+n3)]
      covmat_p[j*n3 +  seq(n2+1,n2+n3) , j*n3 + seq(n2+1+n3*(i-1), n2+n3*(i))] =   covmat_p[seq(n2+1+n3*(i-1), n2+n3*(i)) , seq(n2+1,n2+n3)]
    }
  }
  #chol(covmat_p[-seq(n2), -seq(n2)])
  #image(covmat_p)
  
  # parents current time * parents previous time 
  h = fields::rdist(parents_current_time, parents_previous_times)
  h[h==0] = min(h[h!=0])*.00001
  var_idx = position_in_lower_tri(outer(parents_current_time_var_tag, rep(1, n3)), outer(rep(1, n2), c(parents_previous_times_var_tag)), n_var)
  for(i in seq(2, length(u)))
  {
    covmat_p[seq(n2), seq(n2+(i-2)*n3+1, n2+(i-1)*n3)] = 
      Matern(h, r = matrix(R[i, var_idx], n2), nu_ = matrix(nu_[var_idx], n2)) * 
      sig[i, var_idx]  * rho[var_idx]
    covmat_p[seq(n2+(i-2)*n3+1, n2+(i-1)*n3) , seq(1,n2)] = t(covmat_p[seq(1,n2), seq(n2+(i-2)*n3+1, n2+(i-1)*n3)])
  }
  # checking pdness
  #chol(covmat_p)
  
  ## checking pdness
  #chol(
  #rbind(
  #cbind(covmat_c, covmat_cp),
  #cbind(t(covmat_cp), covmat_p)
  #)
  #)
  #image(covmat_p)
  
  solve_covmat_p = chol2inv(chol(covmat_p))
  salt = covmat_cp %*% solve_covmat_p 
  
}


locs = cbind(runif(100), runif(100))


GMA = function(
    locs, var_tag,
    R, sig,
    rho, nu, 
    n_var
)
{
  #getting cross-variable smoothness from margnial smoothness parameters
  nu_ = outer(nu, nu, "+"); nu_ = nu_/2; nu_ = nu_[lower.tri(nu_, diag = T)]
  # size of spatial sets of interest
  n = nrow(locs)
  #covriance matrix
  covmat = matrix(0, n*length(u), n*length(u))
  # var combinations and distances between unique parent pairs
  var_idx = position_in_lower_tri_cross_vec(var_tag, n_var, diag = T)
  h = as.matrix(dist(locs, diag = T));h = h[lower.tri(h, diag = T)]; 
  h[h==0] = min(h[h!=0])*.00001
  idx = lower_tri_idx(n, diag = T);idx[,1] = idx[,1]
  for(i in seq(length(u))){
    # index of position in the covariance matrix
    idx_ = idx; idx_[,1] = idx_[,1]+(i-1)*n
    covmat[idx_] = 
      Matern(h, r = R[i,var_idx], nu_ = nu_[var_idx]) * 
      sig[i,var_idx] * rho[var_idx]
    # filling symmetric of submatrix
    idx__ = idx[,c(2, 1)]; idx__[,1] = idx__[,1] + n*(i-1)
    covmat[idx__] = covmat[idx_]
    # filling symmetric in block Toeplitz matrix
    if(i>1)covmat[seq(1,n), seq(1+n*(i-1), n*(i))] = covmat[seq(1+n*(i-1), n*(i)) , seq(1,n)] 
    if(i!=(length(u)))for(j in seq(length(u)-i))
    {
      covmat[j*n +  seq(1+n*(i-1), n*(i)) , j*n +  seq(1,n)]                      =   covmat[seq(n*(i-1) + 1, n*(i)) , seq(1,n)]
      if(i!=1)covmat[j*n +  seq(1,n)              , j*n +  seq(1+n*(i-1), n*(i))] =   covmat[seq(n*(i-1) + 1, n*(i)) , seq(1,n)]
    }
  }
  covmat
}




GMA_compressed = function(
    locs, var_tag,
    R, sig,
    rho, nu, 
    n_var
)
{
  #getting cross-variable smoothness from margnial smoothness parameters
  nu_ = outer(nu, nu, "+"); nu_ = nu_/2; nu_ = nu_[lower.tri(nu_, diag = T)]
  # size of spatial sets of interest
  n = nrow(locs)
  #covriance matrix
  covmat = matrix(0, length(u)*n, n)
  # var combinations and distances between unique parent pairs
  var_idx = position_in_lower_tri_cross_vec(var_tag, n_var, diag = T)
  h = as.matrix(dist(locs, diag = T));h = h[lower.tri(h, diag = T)]; 
  h[h==0] = min(h[h!=0])*.00001
  idx = lower_tri_idx(n, diag = T);idx[,1] = idx[,1]
  for(i in seq(length(u))){
    # index of position in the covariance matrix
    idx_ = idx; idx_[,1] = idx_[,1]+(i-1)*n
    covmat[idx_] = 
      Matern(h, r = R[i,var_idx], nu_ = nu_[var_idx]) * 
      sig[i,var_idx] * rho[var_idx]
  }
  covmat
}

expand_covmat = function(compressed_covmat)
{
  n = ncol(compressed_covmat)
  k = nrow(compressed_covmat)/n
  #covriance matrix
  covmat = matrix(0, n*k, n*k)
  covmat[,seq(n)]=compressed_covmat
  idx = lower_tri_idx(n, diag = T);idx[,1] = idx[,1]
  for(i in seq(k)){
    # index of position in the covariance matrix
    idx_ = idx; idx_[,1] = idx_[,1]+(i-1)*n
    # filling symmetric of submatrix
    idx__ = idx[,c(2, 1)]; idx__[,1] = idx__[,1] + n*(i-1)
    covmat[idx__] = covmat[idx_]
    # filling symmetric in block Toeplitz matrix
    if(i>1)covmat[seq(1,n), seq(1+n*(i-1), n*(i))] = covmat[seq(1+n*(i-1), n*(i)) , seq(1,n)] 
    if(i!=(k))for(j in seq(k-i))
    {
      covmat[j*n +  seq(1+n*(i-1), n*(i)) , j*n +  seq(1,n)]                      =   covmat[seq(n*(i-1) + 1, n*(i)) , seq(1,n)]
      if(i!=1)covmat[j*n +  seq(1,n)              , j*n +  seq(1+n*(i-1), n*(i))] =   covmat[seq(n*(i-1) + 1, n*(i)) , seq(1,n)]
    }
  }
  covmat
}




conditional_coeffs = function(precision_bit)
{
  n1 = nrow(precision_bit)
  n2 = ncol(precision_bit)-n1
  #chol_children_precision = chol(precision_bit[,seq(n1)])[seq(n1, 1), seq(n1, 1)]
  #cbind(solve(chol_children_precision) %*% precision_bit[,-seq(n1), drop = F][seq(n1, 1), seq(n2, 1), drop = F], chol_children_precision)
  chol_children_precision = chol(precision_bit[,seq(n1)])
  list("children" = chol_children_precision, "parents" = solve(chol_children_precision) %*% precision_bit[,-seq(n1), drop = F])
}

M = GpGp::exponential_isotropic(c(1, 1, 0), cbind(seq(10), seq(10)))
NNarray = GpGp::find_ordered_nn(cbind(seq(10), seq(10)), 10)
R = GpGp::vecchia_Linv(covparms = c(1, 1, 0), covfun_name = "exponential_isotropic", locs = cbind(seq(10), seq(10)), NNarray = NNarray)
precision_bit = solve(M)[c(1, 2),,drop =F]
conditional_coeffs(precision_bit)

n_children = 10

GMA_vecchia_row = function(
    locs, var_tag,
    R_and_variation, sig_and_variation,
    rho, nu, 
    n_var, 
    n_children, 
    compute_derivative = T
){
  if(n_children>nrow(locs))stop("GMA_Vecchia: there are more children than spatial locations")
  compressed_covmat = GMA_compressed(locs = locs, var_tag = var_tag, R = R_and_variation$R, sig = sig_and_variation$sig, 
               rho = rho, nu = nu, n_var = n_var)
  covmat = expand_covmat(compressed_covmat)
  precision = chol2inv(chol(covmat))
  precision_bit = precision[seq(n_children),,drop = F]
  res = list()
  res$vecchia_row = conditional_coeffs(precision_bit)
  if(compute_derivative)
  {
    nu_variation = array(dim = c(1, length(nu), dim(sig_and_variation[[2]])[3]))
    dimnames(nu_variation)= list(c(""), paste("nu", seq(length(nu))), dimnames(sig_and_variation[[2]])[[3]])
    nu_variation[]= nu; nu_variation[,,seq(dim(nu_variation)[3]-length(nu)+1, dim(nu_variation)[3])] = nu_variation[,,seq(dim(nu_variation)[3]-length(nu)+1, dim(nu_variation)[3])]+diag(.0001, length(nu))
    res$vecchia_row_derivatives = list()
    
    t1 = Sys.time()
    for(i in seq(dim(sig_and_variation[[2]])[3]))
    {
      name = dimnames(sig_and_variation[[2]])[[3]][i]
      print(name)
      # variation of covariance matrix with one parameter
      covmat_variation = expand_covmat(
        GMA_compressed(locs = locs, var_tag = var_tag, R = R_and_variation$R_variation[,,i], sig = sig_and_variation$sig_variation[,,i], 
                                         rho = rho, nu = nu_variation[,,i], n_var = n_var)
        - compressed_covmat
      )
        
      #covmat_variation = as(covmat_variation, "sparseMatrix")
      
      vecchia_row_variation = 
        conditional_coeffs(
          precision_bit +
            as.matrix((precision_bit %*% covmat_variation) %*% precision) # variation of precision following inverse matrix differentiation formula
        )
      res$vecchia_row_derivatives[[name]]=list()
      res$vecchia_row_derivatives[[name]][["children"]]= 10000*(vecchia_row_variation$children-res$vecchia_row$children)
      res$vecchia_row_derivatives[[name]][["parents"]]= 10000*(vecchia_row_variation$parents-res$vecchia_row$parents)
    }
    Sys.time()-t1
    
  }
 
}
image(as.matrix(covmat_variation!=0))


chol(covmat, T)
dim(((precision_bit %*% covmat_variation) %*% precision) [,1:100])
(((precision_bit %*% covmat_variation) %*% precision) [,1:100])
