# seems good except problem with nu when d approx 0


for(i in seq(16))
{
  #################
  # Generate data #
  #################
  print(i)
  # generate test data
  set.seed(i)
  n_var = 10
  n = 40
  
  locs = 10*cbind(runif(n), runif(n))
  var_tag = 1+floor(n_var*runif(n))
  #locs = cbind(runif(30), runif(30))
  #var_tag = 1+rbinom(30, n_var, 1/n_var)
  
  rho = GpGp::exponential_isotropic(c(1, 1, 0), .1*matrix(rnorm(2*n_var), n_var))
  rho_vec = rho[lower.tri(rho, diag = F)];  remove(rho)
  
  
  
  a2_vec = runif(n_var)
  nu_vec = .5 + 2*runif(n_var)
  alpha = .0001 * runif(1)
  a = runif(1) 
  b = runif(1) 
  cc = .1 * runif(1) 
  lambda = runif(1) 
  delta = runif(1) 
  r = .1 * runif(1) 
  A_vec = runif(n_var)
  u = seq(0, 5)
  
  
  ###############################################################################
  # functions to get components of GMA and their variations wrt some parameters #
  ###############################################################################
  
  get_multiplier = function(
    a, b, cc, delta, lambda, r, 
    A_vec, nu_vec, a2_vec, 
    u
  ){
    A_ = outer(A_vec, A_vec, "*") ; A_ = A_[lower.tri(A_, T)]
    a_ = outer(a2_vec, a2_vec, "+") ; a_ =a_/2; a_ = a_[lower.tri(a_, T)]; a_ = sqrt(a_)
    nu_ = outer(nu_vec, nu_vec, "+") ; nu_ = nu_/2; nu_ = nu_[lower.tri(nu_, T)]
    tau_ = a2_vec^(nu_vec)*(1-A_vec^2)^2 / gamma(nu_vec) # n_var 18 τii = σi2 a2νii (1 − Ai ) (here sigma is part of rho, added later)
    tau_ = outer(sqrt(tau_), sqrt(tau_)); tau_ = tau_[lower.tri(tau_, T)]
    res=  
   (
     outer(
       rep(1, length(u)), 
       (2^(1-nu_)) * 
       tau_
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
      A_vec, nu_vec, a2_vec, 
      u
  )
  {
    res= list()
    res$multiplier = get_multiplier(
      a, b, cc, delta, lambda, r, 
      A_vec, nu_vec, a2_vec, 
      u
    )
    res$d_A_vec   = array(0, dim = c(dim(res[[1]]), length(A_vec)))
    res$d_a2_vec   = array(0, dim = c(dim(res[[1]]), length(A_vec)))
    res$d_nu_vec  = array(0, dim = c(dim(res[[1]]), length(A_vec)))
    for(i in seq(A_vec))
    {
      a2_vec_  = a2_vec;  a2_vec_[i]  = a2_vec_[i]+.0001;  res$d_a2_vec[,,i]  = get_multiplier(a, b, cc, delta, lambda, r, A_vec,  nu_vec,  a2_vec_, u)
      A_vec_   = A_vec;   A_vec_[i]   = A_vec_[i]+.0001;   res$d_A_vec[,,i]   = get_multiplier(a, b, cc, delta, lambda, r, A_vec_, nu_vec,  a2_vec,  u)
      nu_vec_  = nu_vec;  nu_vec_[i]  = nu_vec_[i]+.0001;  res$d_nu_vec[,,i]  = get_multiplier(a, b, cc, delta, lambda, r, A_vec,  nu_vec_, a2_vec,  u)
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
  
  i_vec = 1 + floor(n_var*runif(10))
  position_in_lower_tri_cross_vec(i_vec, n_var)
  length(position_in_lower_tri_cross_vec(i_vec, n_var))
  
  
  
  
  
  i_vec = 1 + floor(n_var*runif(10))
  position_in_lower_tri_cross_vec(i_vec, n_var)
  length(position_in_lower_tri_cross_vec(i_vec,n_var))
  
  
  put_ones_in_rho_vec = function(rho_vec)
  {
    nn = .5 + sqrt(.25 + 2*length(rho_vec))
    M = matrix(1, nn, nn)
    M[lower.tri(M)]=rho_vec
    M[lower.tri(M, T)]
  }
  
  ######################
  # getting GMA covmat #
  ######################
  
  
  
  multiplier = get_multiplier(
    a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
    r = r, A_vec = A_vec, nu_vec = nu_vec, a2_vec = a2_vec, u = u
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
    h[h==0] = min(h[h!=0])*.0000001
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
      n_var
    )
  
  #image(expand_covmat(covmat_compressed))
  try(chol(expand_covmat(covmat_compressed)))
  #print(min(eigen(expand_covmat(covmat_compressed))$values))
  print(diag(expand_covmat(covmat_compressed)))
}

########################
# derivative of covmat # 
########################


#rho = GpGp::matern15_isotropic(c(1, 1, 0), matrix(rnorm(2*n_var),n_var))
#rho = diag(1, n_var, n_var)
#rho = matrix(1, n_var, n_var)
#rho_vec = rho[lower.tri(rho, F)];remove(rho)


get_nu_and_variations = function(nu_vec)
{
  res = list()
  nu_ = outer(nu_vec, nu_vec, "+"); nu_= nu_/2; nu_ = nu_[lower.tri(nu_, diag = T)]
  res$nu_ = nu_
  res$nu_variations = matrix(0, length(nu_), length(nu_vec))
  for(j in seq(length(nu_vec)))
  {
    perturb = rep(0, length(nu_vec)); perturb[j]= perturb[j]+.0001
    nu__ = outer(nu_vec + perturb, nu_vec + perturb, "+"); nu__ = nu__/2; nu__ = nu__[lower.tri(nu__, diag = T)]
    res$nu_variations[,j]= nu__
  }
  res
}

nu_and_variations = get_nu_and_variations(nu_vec)




get_rho_and_variations  = function(rho_vec)
{
  res = list()
  res$rho = put_ones_in_rho_vec(rho_vec)
  res$rho_variations = matrix(0, length(res$rho), length(rho_vec))
  for(j in seq(length(rho_vec)))
  {
    perturb = rep(0, length(rho_vec)); perturb[j]= perturb[j]+.0001
    rho_vec_ = rho_vec + perturb
    res$rho_variations[,j]= put_ones_in_rho_vec(rho_vec_)  
  }
  res
}

rho_and_variations = get_rho_and_variations(rho_vec)



multiplier_and_variations = get_multiplier_and_variations(
  a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
  r = r, A_vec = A_vec, nu_vec = nu_vec, a2_vec = a2_vec, u = u
)

effective_range_and_variations = get_effective_range_and_variations(
  a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
  r = r, A_vec = A_vec, a2_vec = a2_vec, u = u
)

#res$covmat_compressed = 
#  GMA_compressed(
#    locs, var_tag,
#    multiplier, 
#    effective_range,
#    nu_vec, 
#    n_var
#  )



GMA_and_derivatives_compressed = function(
    locs, var_tag,
    multiplier_and_variations, 
    effective_range_and_variations,
    nu_and_variations, rho_vec,
    n_var
)
{
  # ading diagonal terms in rho
  rho_vec_ = put_ones_in_rho_vec(rho_vec)
  # size of spatial sets of interest
  n = nrow(locs)
  # var combinations
  var_idx = position_in_lower_tri_cross_vec(var_tag, n_var, diag = T)
  # distance between locs pairs
  h = as.matrix(dist(locs, diag = T));h = h[lower.tri(h, diag = T)]; 
  h[h==0] = min(h[h!=0])*.0000001
  # indices of pairs in covariance matrix
  idx = lower_tri_idx(n, diag = T);idx[,1] = idx[,1]
  # creating matrices
  res = list()
  res$covmat_without_rho = matrix(0, nrow(multiplier)*n, n)
  res$covmat = matrix(0, nrow(multiplier)*n, n)
  res$d_A_vec   = array(0, dim = c(dim(res[[1]]), length(n_var)))
  res$d_a2_vec  = array(0, dim = c(dim(res[[1]]), length(n_var)))
  res$d_nu_vec  = array(0, dim = c(dim(res[[1]]), length(n_var)))
  res$d_rho_vec = array(0, dim = c(dim(res[[1]]), length(rho_vec)))
  
  # looping over time lags
  for(i in seq(nrow(multiplier_and_variations$multiplier))){
    # index of position in the covariance matrix
    idx_ = idx; idx_[,1] = idx_[,1]+(i-1)*n
    # covariance without correlation
    res$covmat_without_rho[idx_] = 
      Matern(h, r = effective_range_and_variations$effective_range[i,var_idx], nu_ = nu_and_variations$nu_[var_idx]) * 
      multiplier_and_variations$multiplier[i, var_idx]
    # covariance
    res$covmat[idx_] = covmat_without_rho[idx_] * rho_vec_[var_idx]
  }
}


t1 = Sys.time()
tatato = lapply(seq(3000), function(i)  covmat_compressed = 
                  GMA_compressed(
                    locs, var_tag,
                    multiplier, 
                    effective_range,
                    nu_vec, 
                    n_var
                  ))
Sys.time()-t1

covmat_compressed = 
  GMA_compressed(
    locs, var_tag,
    multiplier, 
    effective_range,
    nu_vec, 
    n_var
  )
t1 = Sys.time()
tatato = lapply(seq(3000), function(i)  
{
  covmat_compressed = 
    GMA_compressed(
      locs, var_tag,
      multiplier, 
      effective_range,
      nu_vec, 
      n_var
    )
  #for(i in seq(1))expand_covmat(covmat_compressed)
  for(i in seq(1))M = matrix(0, 249, 240)
}
  
  )
Sys.time()-t1





for(i in seq(nrow(multiplier))){
  for(j in seq(dim(multiplier_and_variations$d_A_vec)[3]))
  {
    # index of position in the covariance matrix
    idx_ = idx; idx_[,1] = idx_[,1]+(i-1)*n
    # indices of coefficients that are affected by parameter change
    non_null_idx = which(multiplier_and_variations$d_A_vec[i, var_idx, j]!=multiplier_and_variations$multiplier[i, var_idx])
    res$d_A_vec[,,j][idx_][non_null_idx] = 
      (
        Matern(h[non_null_idx], r = effective_range_and_variations$d_A_vec[i, var_idx, j][non_null_idx], nu_ = nu_and_variations$nu_[var_idx][non_null_idx]) * 
          multiplier_and_variations$d_A_vec[i, var_idx, j][non_null_idx] - 
          res$covmat_compressed[idx_][non_null_idx]
      ) * 10000
    res$d_a2_vec[,,j][idx_][non_null_idx] = 
      (
        Matern(h[non_null_idx], r = effective_range_and_variations$d_a2_vec[i, var_idx, j][non_null_idx], nu_ = nu_and_variations$nu_[var_idx][non_null_idx]) * 
          multiplier_and_variations$d_a2_vec[i, var_idx, j][non_null_idx] - 
          res$covmat_compressed[idx_][non_null_idx]
      ) * 10000
    res$d_nu_vec[,,j][idx_][non_null_idx] = 
      (
        Matern(h[non_null_idx], r = effective_range_and_variations$effective_range[i, var_idx][non_null_idx], nu_ = nu_and_variations$nu_variations[var_idx, j][non_null_idx]) * 
          multiplier_and_variations$d_nu_vec[i, var_idx, j][non_null_idx] - 
          res$covmat_compressed[idx_][non_null_idx]
      ) * 10000
    }
}



multiplier_ = get_multiplier(
  a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
  r = r, A_vec = A_vec, nu_vec = nu_vec + c(.0001, rep(0, 9)), a2_vec = a2_vec,u = u
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
  effective_range,
  nu_vec_, 
  n_var
)-
res$covmat_compressed 
  )
)







  
  
  
  


