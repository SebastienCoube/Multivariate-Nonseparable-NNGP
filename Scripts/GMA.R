
for(i in seq(10))
{
  print(i)
# generate test data
set.seed(i)
p = 10
n = 40

rho = GpGp::matern15_isotropic(c(1, 1, 0), matrix(rnorm(2*p),p))
rho = diag(1, p, p)
rho = matrix(1, p, p)
rho_vec = rho[lower.tri(rho, T)];remove(rho)

locs = 10*cbind(runif(n), runif(n))
var_tag = 1+floor(p*runif(n))
#locs = cbind(runif(30), runif(30))
#var_tag = 1+rbinom(30, p, 1/p)

a2_vec = runif(p)
nu_vec = .5 + 2*runif(p)
alpha = .0001 * runif(1)
a = runif(1) 
b = runif(1) 
cc = .1 * runif(1) 
lambda = runif(1) 
delta = runif(1) 
r = .1 * runif(1) 
A_vec = runif(p)
u = seq(0, 10)


###############################################################################
# functions to get components of GMA and their variations wrt some parameters #
###############################################################################

get_multiplier = function(
  a, b, cc, delta, lambda, r, 
  A_vec, nu_vec, a2_vec, rho_vec, 
  u
){
  A_ = outer(A_vec, A_vec, "*") ; A_ = A_[lower.tri(A_, T)]
  a_ = outer(a2_vec, a2_vec, "+") ; a_ =a_/2; a_ = a_[lower.tri(a_, T)]; a_ = sqrt(a_)
  nu_ = outer(nu_vec, nu_vec, "+") ; nu_ = nu_/2; nu_ = nu_[lower.tri(nu_, T)]
  tau_ = a2_vec^(nu_vec)*(1-A_vec^2)^2 / gamma(nu_vec) # p 18 τii = σi2 a2νii (1 − Ai ) (here sigma is part of rho)
  tau_ = outer(sqrt(tau_), sqrt(tau_)); tau_ = tau_[lower.tri(tau_, T)]
  res=  
 (
   outer(
     rep(1, length(u)), 
     (2^(1-nu_)) * 
     tau_ * 
      rho_vec
     * #gamma(nu_)*
       a_^(-2*nu_)
   )
 ) / (
      (
        outer((1 + (cc*u)^(2*a)     )^  delta,  rep(1, length(A_))) - 
        outer((1 + (r *u)^(2*lambda))^(-delta), A_                )
      )
      * 
        (
          outer((1 + (cc*u)^(2*a)     )^  b,  rep(1, length(A_))) - 
          outer((1 + (r *u)^(2*lambda))^(-b), A_                )
        )
    )
  row.names(res) = paste("u=", u, sep="")
  M = matrix(0, length(A_vec),length(A_vec))
  colnames(res) = paste("var", row(M)[lower.tri(M, T)], col(M)[lower.tri(M, T)] , sep="_")
  res
}

get_multiplier_and_variations = function(
    a, b, cc, delta, lambda, r, 
    A_vec, nu_vec, a2_vec, rho_vec, 
    u
)
{
  res= list()
  res$multiplier = get_multiplier(
    a, b, cc, delta, lambda, r, 
    A_vec, nu_vec, a2_vec, rho_vec, 
    u
  )
  res$d_A_vec   = array(0, dim = c(dim(res[[1]]), length(A_vec)))
  res$d_a2_vec   = array(0, dim = c(dim(res[[1]]), length(A_vec)))
  res$d_nu_vec  = array(0, dim = c(dim(res[[1]]), length(A_vec)))
  #res$d_rho_vec = array(0, dim = c(dim(res[[1]]), length(A_vec)))
  for(i in seq(A_vec))
  {
    a2_vec_  = a2_vec;  a2_vec_[i]  = a2_vec_[i]+.0001;  res$d_a2_vec[,,i]  = get_multiplier(a, b, cc, delta, lambda, r, A_vec,  nu_vec,  a2_vec_, rho_vec, u)
    A_vec_   = A_vec;   A_vec_[i]   = A_vec_[i]+.0001;   res$d_A_vec[,,i]   = get_multiplier(a, b, cc, delta, lambda, r, A_vec_, nu_vec,  a2_vec,  rho_vec, u)
    nu_vec_  = nu_vec;  nu_vec_[i]  = nu_vec_[i]+.0001;  res$d_nu_vec[,,i]  = get_multiplier(a, b, cc, delta, lambda, r, A_vec,  nu_vec_, a2_vec,  rho_vec, u)
    #rho_vec_ = rho_vec; rho_vec_[i] = rho_vec_[i]+.0001; res$d_rho_vec[,,i] = (get_multiplier(a, b, cc, delta, lambda, r, A_vec,  nu_vec,  a2_vec,  rho_vec_, u)-res[[1]])*10000
  }
  res
}

get_effective_range = function(
  a, b, cc, delta, lambda, r, 
  A_vec, a2_vec, 
  u
)
{
  A_ = outer(A_vec, A_vec, "*") ; A_ = A_[lower.tri(A_, T)]
  a_ = outer(a2_vec, a2_vec, "+") ; a_ =a_/2; a_ = a_[lower.tri(a_, T)]; a_ = sqrt(a_)
  outer(
    rep(1, length(u)), 
    a_
    )/
    (
      outer((1 + (cc*u)^(2*a)     )^  b,  rep(1, length(A_))) - 
      outer((1 + (r *u)^(2*lambda))^(-b), A_                )
    )^.5
}



get_effective_range_and_variations = function(
    a, b, cc, delta, lambda, r, 
    A_vec, nu_vec, a2_vec,  
    u
)
{
  res= list()
  res$effective_range = get_effective_range(
    a, b, cc, delta, lambda, r, 
    A_vec, a2_vec,  
    u
  )
  res$d_A_vec   = array(0, dim = c(dim(res[[1]]), length(A_vec)))
  res$d_a2_vec   = array(0, dim = c(dim(res[[1]]), length(A_vec)))
  for(i in seq(A_vec))
  {
    a2_vec_  = a2_vec;  a2_vec_[i]  = a2_vec_[i]+.0001;  res$d_a2_vec[,,i]  = get_effective_range(a, b, cc, delta, lambda, r, A_vec,  a2_vec_, u)
    A_vec_   = A_vec;   A_vec_[i]   = A_vec_[i]+.0001;   res$d_A_vec[,,i]   = get_effective_range(a, b, cc, delta, lambda, r, A_vec_, a2_vec,  u)
  }
  res
}

Matern = function(h, r, nu_)(r * h)^nu_ * besselK(r * h, nu_)# / gamma(nu_)

##################################################
# Some functions to work with symmetric matrices #
##################################################

# given integer n, gives 2-column matrix of lower triangular indices 
# in the square matrix n*n WITHOUT diagonal as in dist(,diag = F)


# ex :              index | i j
# n= 4 -> 1 . . .    1    | 1 1
#         2 5 . .    2    | 2 1
#         3 6 8 .    3    | 3 1
#         4 7 9 10   4    | 4 1
#                    5    | 2 2
#                    6    | 3 2 
#                    7    | 4 2 
#                    8    | 3 3 
#                    9    | 4 3
#                    10   | 4 4 
#                      
    
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

# ex: i = 4, j = 3, n = 4
#             j
#         1 . . . 
#         2 5 . . 
#         3 6 8 . 
#       i 4 7(9)10
#     result = 9

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





n_var = 3
i_vec = 1 + floor(3*runif(10))
position_in_lower_tri_cross_vec(i_vec, 3)
length(position_in_lower_tri_cross_vec(i_vec, 3))

######################
# getting GMA covmat #
######################



multiplier = get_multiplier(
  a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
  r = r, A_vec = A_vec, nu_vec = nu_vec, a2_vec = a2_vec, rho_vec = rho_vec, u = u
  )
effective_range = get_effective_range(
  a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
  r = r, A_vec = A_vec, a2_vec = a2_vec, u = u
  )


GMA_compressed = function(
    locs, var_tag,
    multiplier, 
    effective_range,
    nu_vec, 
    n_var
)
{
  #getting cross-variable smoothness from margnial smoothness parameters
  nu_ = outer(nu_vec, nu_vec, "+"); nu_ = nu_/2; nu_ = nu_[lower.tri(nu_, diag = T)]
  # size of spatial sets of interest
  n = nrow(locs)
  #covriance matrix
  covmat = matrix(0, nrow(multiplier)*n, n)
  # var combinations and distances between unique parent pairs
  var_idx = position_in_lower_tri_cross_vec(var_tag, n_var, diag = T)
  h = as.matrix(dist(locs, diag = T));h = h[lower.tri(h, diag = T)]; 
  h[h==0] = min(h[h!=0])*.00001
  idx = lower_tri_idx(n, diag = T);idx[,1] = idx[,1]
  for(i in seq(nrow(multiplier))){
    # index of position in the covariance matrix
    idx_ = idx; idx_[,1] = idx_[,1]+(i-1)*n
    covmat[idx_] = 
      Matern(h, r = effective_range[i,var_idx], nu_ = nu_[var_idx]) * 
      multiplier[i, var_idx]
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

covmat_compressed = 
  GMA_compressed(
    locs, var_tag,
    multiplier, 
    effective_range,
    nu_vec, 
    p
  )

#image(expand_covmat(covmat_compressed))
try(chol(expand_covmat(covmat_compressed)))
#print(min(eigen(expand_covmat(covmat_compressed))$values))
#print(diag(expand_covmat(covmat_compressed)))
}


multiplier_and_variations = get_multiplier_and_variations(
  a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
  r = r, A_vec = A_vec, nu_vec = nu_vec, a2_vec = a2_vec, rho_vec = rho_vec, u = u
)



effective_range_and_variations = get_effective_range_and_variations(
  a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
  r = r, A_vec = A_vec, a2_vec = a2_vec, u = u
)

res = list()
res$covmat_compressed = 
  GMA_compressed(
    locs, var_tag,
    multiplier, 
    effective_range,
    nu_vec, 
    p
  )

n_var = p

res$d_A_vec   = array(0, dim = c(dim(res[[1]]), length(A_vec)))
res$d_a2_vec   = array(0, dim = c(dim(res[[1]]), length(A_vec)))
res$d_nu_vec  = array(0, dim = c(dim(res[[1]]), length(A_vec)))
res$d_rho_vec  = array(0, dim = c(dim(res[[1]]), length(rho_vec) - length(A_vec)))


#getting cross-variable smoothness from margnial smoothness parameters
nu_ = outer(nu_vec, nu_vec, "+"); nu_ = nu_/2; nu_ = nu_[lower.tri(nu_, diag = T)]
# size of spatial sets of interest
n = nrow(locs)
# var combinations and distances between unique parent pairs
var_idx = position_in_lower_tri_cross_vec(var_tag, n_var, diag = T)
h = as.matrix(dist(locs, diag = T));h = h[lower.tri(h, diag = T)]; 
h[h==0] = min(h[h!=0])*.00001
idx = lower_tri_idx(n, diag = T);idx[,1] = idx[,1]

for(j in seq(dim(multiplier_and_variations$d_A_vec)[3]))
{
  for(i in seq(nrow(multiplier))){
    # index of position in the covariance matrix
    idx_ = idx; idx_[,1] = idx_[,1]+(i-1)*n
    # indices of coefficients that are affected by parameter change
    non_null_idx = which(multiplier_and_variations$d_A_vec[i, var_idx, j]!=multiplier_and_variations$multiplier[i, var_idx])
    res$d_A_vec[,,j][idx_][non_null_idx] = 
      (
        Matern(h[non_null_idx], r = effective_range_and_variations$d_A_vec[i, var_idx, j][non_null_idx], nu_ = nu_[var_idx][non_null_idx]) * 
          multiplier_and_variations$d_A_vec[i, var_idx, j][non_null_idx] - 
          res$covmat_compressed[idx_][non_null_idx]
      ) * 10000
    res$d_a2_vec[,,j][idx_][non_null_idx] = 
      (
        Matern(h[non_null_idx], r = effective_range_and_variations$d_a2_vec[i, var_idx, j][non_null_idx], nu_ = nu_[var_idx][non_null_idx]) * 
          multiplier_and_variations$d_a2_vec[i, var_idx, j][non_null_idx] - 
          res$covmat_compressed[idx_][non_null_idx]
      ) * 10000
    perturb = rep(0, length(nu_vec)); perturb[j]= perturb[j]+.0001
    nu__ = outer(nu_vec + preturb, nu_vec + preturb, "+"); nu__ = nu__/2; nu__ = nu__[lower.tri(nu_, diag = T)]
    res$d_nu_vec[,,j][idx_][non_null_idx] = 
      (
        Matern(h[non_null_idx], r = effective_range_and_variations$effective_range[i, var_idx][non_null_idx], nu_ = nu__[var_idx][non_null_idx]) * 
          multiplier_and_variations$d_nu_vec[i, var_idx, j][non_null_idx] - 
          res$covmat_compressed[idx_][non_null_idx]
      ) * 10000
    }
}



multiplier_ = get_multiplier(
  a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
  r = r, A_vec = A_vec, nu_vec = nu_vec + c(.0001, rep(0, 9)), a2_vec = a2_vec, rho_vec = rho_vec, u = u
)

effective_range_ = get_effective_range(
  a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
  r = r, A_vec = A_vec, a2_vec = a2_vec #+ c(.0001, rep(0, 9))
  , u = u
)

nu_vec_ = nu_vec
nu_vec_[1]=nu_vec_[1]+.0001

plot(
  
res$d_nu_vec[,,1],
  
10000*(
GMA_compressed(
  locs, var_tag,
  multiplier_, 
  effective_range_,
  nu_vec_, 
  p
)-
res$covmat_compressed 
  )
)
  
  
  
  


