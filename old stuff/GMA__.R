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
  n_locs = 40
  
  locs = 10*cbind(runif(n_locs), runif(n_locs))
  var_tag = 1+floor(n_var*runif(n_locs))
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
  
  # given integer n_locs, gives 2-column matrix of lower triangular indices 
  # in the square matrix n_locs*n_locs WITHOUT diagonal as in dist(,diag = F)
  
  
  # ex :              index | i j
  # n_locs= 4 -> 1 . . .    1    | 1 1
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
  
  lower_tri_idx = function(n_locs, diag = F)
  {
    
    if(diag == F)
    {
      return(
        cbind(
          rev(abs(sequence(seq.int(n_locs - 1)) - n_locs) + 1),
          rep.int(seq.int(n_locs - 1), rev(seq.int(n_locs - 1)))
        )
      )
    }
    if(diag == T)
    {
      return(
        cbind(
          rev(abs(sequence(seq.int(n_locs)) - n_locs) + 1),
          rep.int(seq.int(n_locs), rev(seq.int(n_locs)))
        )
      )
    }
  }
  lower_tri_idx(10)
  
  # given three integers i, j, n_locs, gives the position of coefficient (i, j)
  # in the n_locs*n_locs square matrix in the lower triangular coefficients
  # WITH diagonal as in dist(,diag = T)
  
  # ex: i = 4, j = 3, n_locs = 4
  #             j
  #         1 . . . 
  #         2 5 . . 
  #         3 6 8 . 
  #       i 4 7(9)10
  #     result = 9
  
  position_in_lower_tri = function(i, j, n_locs)
  {  
    a = pmin(i, j); b = pmax(i, j)
    (a-1)*(n_locs-a+1) + (a-1)*a/2 + b-a+1
  }
  #  i = 1;j =3;n_locs=3;position_in_lower_tri(i, j, n_locs)
  #  i = 2;j =1;n_locs=3;position_in_lower_tri(i, j, n_locs)
  #  i = 2;j =3;n_locs=3;position_in_lower_tri(i, j, n_locs)
  #  i = 3;j =3;n_locs=3;position_in_lower_tri(i, j, n_locs)
  
  
  # given a vector of integers i_vec and an integer n_var, 
  # n_locs var being greater than the entries of i_vec, 
  # computes the lower-triangular without diagonal of the crossing matrix of i_vec, i_vec
  # and finds the position of those pairs in the lower triangular part of the 
  # matrix of size n_var * n_var
  position_in_lower_tri_cross_vec = function(i_vec, n_var, diag = F)
  {  
    idx = lower_tri_idx(length(i_vec), diag = diag)
    position_in_lower_tri(i_vec[idx[,1]], i_vec[idx[,2]], n_var)
  }
  
##  i_vec = 1 + floor(n_var*runif(10))
##  position_in_lower_tri_cross_vec(i_vec, n_var)
##  length(position_in_lower_tri_cross_vec(i_vec, n_var))
##  
##  i_vec = 1 + floor(n_var*runif(10))
##  position_in_lower_tri_cross_vec(i_vec, n_var)
##  length(position_in_lower_tri_cross_vec(i_vec,n_var))
  
  
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
    n_var, 
    rho_vec_with_ones
  )
  {
    #getting cross-variable smoothness from margnial smoothness parameters
    nu_ = outer(nu_vec, nu_vec, "+"); nu_ = nu_/2; nu_ = nu_[lower.tri(nu_, diag = T)]
    # size of spatial sets of interest
    n_locs = nrow(locs)
    #covriance matrix
    covmat = matrix(0, nrow(multiplier)*n_locs, n_locs)
    # var combinations and distances between unique parent pairs
    var_idx = position_in_lower_tri_cross_vec(var_tag, n_var, diag = T)
    h = as.matrix(dist(locs, diag = T));h = h[lower.tri(h, diag = T)]; 
    h[h==0] = min(h[h!=0])*.0000001
    covmat_coeffs = 
      Matern(h, r = t(effective_range)[var_idx,], nu_ = nu_[var_idx]) * 
      t(multiplier)[var_idx,] * rho_vec_with_ones[var_idx]
    covmat_coeffs
  }
  
  covmat_coeffs  =   GMA_compressed(
    locs, var_tag,
    multiplier, 
    effective_range,
    nu_vec, 
    n_var, 
    put_ones_in_rho_vec(rho_vec)
  )
  block_lower_tri_idx = lower_tri_idx(n_locs, T)
  expand_covmat_into_blocks = function(covmat_coeffs, n_locs, block_lower_tri_idx)
  {
    #if(!sparse)
      blocks = lapply(seq(ncol(covmat_coeffs)), function(i)
    {
      M = matrix(0, n_locs, n_locs)
      M[block_lower_tri_idx]= covmat_coeffs[, i]
      M = t(M)
      M[block_lower_tri_idx]= covmat_coeffs[, i]
      M
    })
    blocks  
  }
  
  expand_block_toeplitz_covmat = function(covmat_coeffs, n_locs, block_lower_tri_idx)
  {
    blocks = expand_covmat_into_blocks(covmat_coeffs, n_locs, block_lower_tri_idx)
    block_idx_matrix = toeplitz(c(seq(ncol(covmat_coeffs))))
    res=  matrix(0, n_locs*ncol(covmat_coeffs), n_locs*ncol(covmat_coeffs))
    for(i in seq(ncol(block_idx_matrix)))
    {
      for(j in seq(ncol(block_idx_matrix))){
        res[seq(((i-1)*n_locs+1), (i)*n_locs), seq(((j-1)*n_locs+1),(j)*n_locs)]=blocks[[block_idx_matrix[i,j]]]
      }
    }
    res
  }

  
  Rcpp::sourceCpp("multiply.cpp")
  #####  
  #####  multiply_vector_block_toeplitz = function(v, covmat_coeffs, n_locs, idx_mat)
  #####  {
  #####    block_idx_matrix = toeplitz(c(seq(ncol(covmat_coeffs))))
  #####    res=  rep(0, n_locs*ncol(covmat_coeffs))
  #####    for(i in seq(ncol(block_idx_matrix)))
  #####    {
  #####      for(j in seq(ncol(block_idx_matrix)))
  #####      {
  #####        res[(1+(i-1)*n_locs):(i*n_locs)] = res[(1+(i-1)*n_locs):(i*n_locs)] + Multiply_vector_compressed_symmat(
  #####          v = v[(1+(j-1)*n_locs):(j*n_locs)], compressed_symmat = covmat_coeffs[,block_idx_matrix[j, i]], 
  #####          idx_mat = idx_mat, n_locs = n_locs
  #####          )
  #####      }
  #####    }
  #####    res
  #####  }
  #####  
  #####  multiply_vector_block_toeplitz = function(v, covmat_coeffs, n_locs, block_lower_tri_idx)
  #####  {
  #####    # looking for null coefficients in covmat coeffs
  #####    non_zero_idx  =which(covmat_coeffs[,1]!=0)
  #####    # removing null coefficients
  #####    idx = block_lower_tri_idx[non_zero_idx,]
  #####    covmat_coeffs[-non_zero_idx,]=NULL
  #####    block_lower_tri_idx[-non_zero_idx,]=NULL
  #####    
  #####    M = Matrix::sparseMatrix(i = idx[,2], j= idx[,1], x= covmat_coeffs[,1], symmetric = T)
  #####    
  #####    v = matrix(v, n_locs)
  #####    vM = matrix(0, n_locs, ncol(v))
  #####    
  #####    for(i in seq(ncol(covmat_coeffs)))
  #####    {
  #####      M@x = covmat_coeffs[,i]
  #####      toplitz_idx = Matrix::sparseVector(1, i, ncol(covmat_coeffs))
  #####      vM = vM + M %*% v %*% Matrix::toeplitz(toplitz_idx)
  #####    }
  #####    vM
  #####  }
  
  
  multiply_vector_block_toeplitz_sparse = function(v, covmat_coeffs, n_locs, idx_mat)
  {
    k = length(blocks)
    v_ = matrix(v, ncol = k)
    res = matrix(0, n_locs, k)
    for(lag_idx in seq(k))
    {
      toeplitz_idx = rep(0, k)
      toeplitz_idx[lag_idx]=1
      res = res + Multiply_matrix_compressed_symmat(v_, covmat_coeffs[,lag_idx], idx_mat, n_locs) %*% toeplitz(toeplitz_idx)
    }
    res
  }
  
  multiply_vector_block_toeplitz = function(v, blocks, n_locs)
  {
    k = length(blocks)
    v_ = matrix(v, ncol = k)
    res = matrix(0, n_locs, k)
    for(lag_idx in seq(k))
    {
      toeplitz_idx = rep(0, k)
      toeplitz_idx[lag_idx]=1
      res = res + blocks [[lag_idx]] %*% v_ %*% toeplitz(toeplitz_idx)
    }
    res
  }
  
  #Start from here
##  multiply_d_rho = function(blocks, v, n_locs, var_tag, n_var)
##  {
##    res = array(0, c(n_locs, length(blocks), n_var* (n_var-1)/2))
##    v_ = matrix(v, n_locs)
##    for(lag_idx_1 in seq(length(blocks))){
##      for(lag_idx_2 in seq(length(blocks))){
##        v_blocks= v[,lag_idx_2] %*% blocks[lag_idx_1]
##        for(var_idx_1 in seq(nvar-1)){
##          for(var_idx_2 in seq(var_idx_1+1, nvar)){
##            idx = which(var_tag==var_idx_1|var_tag==var_idx_2)
##            
##        }
##      }
##    }
##    }  
    
  



  
  
  v = rnorm(n_locs*length(u))
  blocks = expand_covmat_into_blocks(covmat_coeffs, n_locs, block_lower_tri_idx)
  multiply_vector_block_toeplitz(v = v, blocks, n_locs)- 
    matrix(expand_block_toeplitz_covmat(covmat_coeffs, n_locs, block_lower_tri_idx) %*% v, n_locs)
  
  
  multiply_vector_block_toeplitz(v = v, blocks, n_locs)- 
    matrix(expand_block_toeplitz_covmat(covmat_coeffs, n_locs, block_lower_tri_idx) %*% v, n_locs)
  
  
  multiply_vector_block_toeplitz_sparse(v = v, covmat_coeffs = covmat_coeffs, idx_mat = block_lower_tri_idx, n_locs)- 
    matrix(expand_block_toeplitz_covmat(covmat_coeffs, n_locs, block_lower_tri_idx) %*% v, n_locs)
  
  
  # idea : multiply block by v, the by a var idx matrix
  # subset the result
  t1 = Sys.time()
  tatato = lapply(seq(36 * 1000), function(x)
  {
    var_idx_matrix = matrix(0, n_var, n_locs)
    var_idx_matrix[cbind(var_tag, seq(n_locs))]=1
    res = matrix(0, n_var*(n_var-1)/2, n_locs)
    left_mult = var_idx_matrix %*% blocks[[1]]
    right_mult = blocks[[1]] %*% t(var_idx_matrix)
    #k = 1
    #for(i in seq(n_var-1))
    #{
    #  for(j in seq(i+1, n_var))
    #  {
    #    idx = which(var_tag==j)
    #    res[k,idx] = res[k,idx] + left_mult[i,idx]
    #    idx = which(var_tag==i)
    #    res[k,idx] = res[k,idx] + right_mult[idx, j]
    #    k=k+1
    #  }
    #}
  }
  )
  Sys.time()-t1
  

    
  
  
  
    t1 =Sys.time()
    tatato= lapply(seq(10000), function(i)
    {
      tatat0= multiply_vector_block_toeplitz(v = v, blocks = expand_covmat_into_blocks(covmat_coeffs, n_locs, block_lower_tri_idx), n_locs=n_locs)
      return(0)
    })
    Sys.time()-t1
  
  
  
    
    t1 =Sys.time()
    tatato= lapply(seq(10000), function(i)
    {
      tatat0= multiply_vector_block_toeplitz(v = v, blocks = expand_covmat_into_blocks(covmat_coeffs[seq(8),], n_locs, block_lower_tri_idx[seq(8),]), n_locs=n_locs)
      return(0)
    })
    Sys.time()-t1
    
    t1 =Sys.time()
    tatato= lapply(seq(10000), function(i)
    {
      tatat0= multiply_vector_block_toeplitz_sparse(v = v, covmat_coeffs = covmat_coeffs, n_locs = n_locs, idx_mat = block_lower_tri_idx[seq(8),])
      return(0)
    })
    Sys.time()-t1
    
    
    
  ##  
  ##  
  ##  t1 =Sys.time()
  ##  tatato= lapply(seq(1000), function(i)
  ##  {
  ##    tatat0= v %*% (expand_block_toeplitz_covmat(covmat_coeffs = 
  ##                                                GMA_compressed(
  ##                                                  locs, var_tag,
  ##                                                  multiplier, 
  ##                                                  effective_range,
  ##                                                  nu_vec, 
  ##                                                  n_var
  ##                                                ), n_locs = n_locs, block_lower_tri_idx = block_lower_tri_idx))
  ##    return(0)
  ##  })
  ##  Sys.time()-t1
  ##  
  ##  
  ##  t1 =Sys.time()
  ##  tatato= lapply(seq(1000), function(i)
  ##  {
  ##    tatat0= chol(expand_block_toeplitz_covmat(covmat_coeffs = 
  ##                                                GMA_compressed(
  ##                                                  locs, var_tag,
  ##                                                  multiplier, 
  ##                                                  effective_range,
  ##                                                  nu_vec, 
  ##                                                  n_var
  ##                                                ), n_locs = n_locs, block_lower_tri_idx = block_lower_tri_idx))
  ##    return(0)
  ##  })
  ##  Sys.time()-t1
  
  
  
  
  image(expand_block_toeplitz_covmat(covmat_coeffs, n_locs, block_lower_tri_idx))
  try(chol(expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs, n_locs = n_locs, block_lower_tri_idx = block_lower_tri_idx)))
  #print(min(eigen(expand_covmat(covmat_compressed))$values))
  print(diag(expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs, n_locs = n_locs, block_lower_tri_idx = block_lower_tri_idx)))
}

########################
# derivative of covmat # 
########################




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


multiplier_and_variations = get_multiplier_and_variations(
  a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
  r = r, A_vec = A_vec, nu_vec = nu_vec, a2_vec = a2_vec, u = u
)

effective_range_and_variations = get_effective_range_and_variations(
  a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
  r = r, A_vec = A_vec, a2_vec = a2_vec, u = u
)



rho_vec_with_ones = put_ones_in_rho_vec(rho_vec)

GMA_and_derivatives_compressed = function(
    locs, var_tag,
    multiplier_and_variations, 
    effective_range_and_variations,
    nu_and_variations, rho_vec_with_ones
)
{
  # size of spatial sets of interest
  n_locs = nrow(locs)
  n_lags = nrow(multiplier_and_variations$multiplier)
  n_var = dim(multiplier_and_variations$d_A_vec)[3]
  # var combinations
  var_idx = position_in_lower_tri_cross_vec(var_tag, n_var, diag = T)
  # distance between locs pairs
  h = as.matrix(dist(locs, diag = T));h = h[lower.tri(h, diag = T)]
  h[h==0] = min(h[h!=0])*.0001
  # indices of pairs in covariance matrix
  idx = lower_tri_idx(n_locs, diag = T);idx[,1] = idx[,1]
  # creating matrices
  res = list()
  res$covmat_without_rho = matrix(0, nrow(idx), n_lags)
  res$covmat = matrix(0, nrow(idx), n_lags)
  res$d_A_vec   = array(0, dim = c(dim(res[[1]]), n_var))
  res$d_a2_vec  = array(0, dim = c(dim(res[[1]]), n_var))
  res$d_nu_vec  = array(0, dim = c(dim(res[[1]]), n_var))
  
  # looping over time lags
  for(i_lag in seq(n_lags)){
    # index of position in the covariance matrix
    # covariance without correlation
    res$covmat_without_rho[,i_lag] = 
      Matern(h, r = effective_range_and_variations$effective_range[i_lag,var_idx], nu_ = nu_and_variations$nu_[var_idx]) * 
      multiplier_and_variations$multiplier[i_lag, var_idx]
    # covariance
    res$covmat[,i_lag] = res$covmat_without_rho[,i_lag] * rho_vec_with_ones[var_idx]
    # variations of covariance
    for(i_var in seq(n_var))
    {
      res$d_A_vec[,i_lag, i_var] = 
        Matern(h, r = effective_range_and_variations$d_A_vec[i_lag, var_idx, i_var], nu_ = nu_and_variations$nu_[var_idx]) * 
        multiplier_and_variations$d_A_vec[i_lag, var_idx, i_var] * rho_vec_with_ones[var_idx]
      res$d_a2_vec[,i_lag, i_var] = 
        Matern(h, r = effective_range_and_variations$d_a2_vec[i_lag, var_idx, i_var], nu_ = nu_and_variations$nu_[var_idx]) * 
        multiplier_and_variations$d_a2_vec[i_lag, var_idx, i_var] * rho_vec_with_ones[var_idx]
      res$d_nu_vec[,i_lag, i_var] = 
        Matern(h, r = effective_range_and_variations$effective_range[i_lag,var_idx], nu_ = nu_and_variations$nu_variations[var_idx, i_var]) * 
        multiplier_and_variations$multiplier[i_lag, var_idx] * rho_vec_with_ones[var_idx]
    }
  }
  res
}



t1 = Sys.time()
tatat0 = sapply(seq(100), function(i)GMA_and_derivatives_compressed (
    locs = locs, var_tag = var_tag,
    multiplier_and_variations = multiplier_and_variations, 
    effective_range_and_variations = effective_range_and_variations,
    nu_and_variations = nu_and_variations, rho_vec_with_ones = rho_vec_with_ones 
))
Sys.time()-t1


GMA_compressed = function(
    locs, var_tag,
    multiplier, 
    effective_range,
    nu_vec, 
    n_var, 
    rho_vec
)
{
  #getting cross-variable smoothness from margnial smoothness parameters
  nu_ = outer(nu_vec, nu_vec, "+"); nu_ = nu_/2; nu_ = nu_[lower.tri(nu_, diag = T)]
  # size of spatial sets of interest
  n_locs = nrow(locs)
  #covriance matrix
  covmat = matrix(0, nrow(multiplier)*n_locs, n_locs)
  # var combinations and distances between unique parent pairs
  var_idx = position_in_lower_tri_cross_vec(var_tag, n_var, diag = T)
  h = as.matrix(dist(locs, diag = T));h = h[lower.tri(h, diag = T)]; 
  h[h==0] = min(h[h!=0])*.0000001
  covmat_coeffs = 
    Matern(h, r = t(effective_range)[var_idx,], nu_ = nu_[var_idx]) * 
    t(multiplier)[var_idx,] * rho_vec[var_idx]
  covmat_coeffs
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

t1 = Sys.time()
tatato = lapply(seq(3000), function(i)  covmat_compressed = 
                  GMA_and_derivatives_compressed(
                    locs, var_tag,
                    multiplier_and_variations = multiplier_and_variations, 
                    effective_range_and_variations = effective_range_and_variations,
                    nu_and_variations = nu_and_variations, 
                    rho_vec = rho_vec,
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


for(i in seq(nrow(multiplier))){
  for(j in seq(dim(multiplier_and_variations$d_A_vec)[3]))
  {
    # index of position in the covariance matrix
    idx_ = idx; idx_[,1] = idx_[,1]+(i-1)*n_locs
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













