
#### # flip image function to match matrix
#### my_image = function(m)image(t(m)[,nrow(m):1])
#### 
#### # simulate data
#### set.seed(2)
#### 
#### n_var = 8
#### time_depth = 4
#### 
#### # spatial locations and variable index
#### locs_ = as.matrix(expand.grid(seq(100), seq(100)))/20#cbind(runif(n_loc), runif(n_loc))
#### locs_ = locs_ + rnorm(length(locs_), 0, .0001)
#### locs_ = locs_[sample(seq(nrow(locs_)), size = nrow(locs_), F),]
#### n_loc =nrow(locs_)
#### row.names(locs_) = seq(n_loc)
#### var_tag = ceiling(n_var*runif(n_loc));var_tag = match(var_tag, unique(var_tag))
#### 
#### # covariance params
#### rho = GpGp::exponential_isotropic(c(1, 1, 0), matrix(rnorm(2*n_var), n_var))
#### rho_vec = rho[lower.tri(rho, diag = F)];  remove(rho)
#### a2_vec = rep(100, n_var)#.000001*runif(n_var)
#### nu_vec = .5 + 2*runif(n_var)
#### alpha = .01 * runif(1)
#### a = runif(1) 
#### b = runif(1) 
#### cc = 1*runif(1) 
#### lambda = runif(1) 
#### delta = runif(1) 
#### r = 1 * runif(1) 
#### A_vec = rep(0, n_var)#runif(n_var)
#### u = seq(0, time_depth-1)
#### 
#### # creating the DAG
#### NNarray_same_time = find_ordered_nn_multi(locs_ = locs_, 
####                                           var_tag = var_tag, 
####                                           m_whatever_closest = 10, m_same_var = 5, m_other_vars = 0, 
####                                           lonlat = F)
#### NNarray_pevious_times = find_unordered_nn_multi(locs_ = locs_, 
####                                                 var_tag = var_tag, 
####                                                 m_whatever_closest = 10, m_same_var = 5, m_other_vars = 0, 
####                                                 lonlat = F)
#### DAG = list(
####   "children" = lapply(seq(n_loc), function(x)x), 
####   "parents_same_time" = NNarray_same_time, 
####   "parents_previous_times" = NNarray_pevious_times 
#### )
#### 

# get lower triangular indices in the dag
get_lower_tri_idx_DAG = function(DAG)
{
  lower_tri_idx_DAG = list()
  lower_tri_idx_DAG$same_time = lapply(mapply(c, DAG$parents_same_time, DAG$children, SIMPLIFY = F), function(x)get_lower_tri_idx(length(x), T))
  if(!is.null(DAG$parents_previous_times))lower_tri_idx_DAG$previous_times = lapply(DAG$parents_previous_times, function(x)get_lower_tri_idx(length(x), T))
  lower_tri_idx_DAG
}


# get var idx in lower tri mat
get_var_idx = function(DAG, var_tag)
{
  n_var = max(var_tag)
  var_idx = list()
  var_idx$current_time = lapply(mapply(c,DAG$parents_same_time, DAG$children, SIMPLIFY = F), function(idx)position_in_lower_tri_cross_vec(var_tag[idx], n_var = n_var, diag = T))
  if(!is.null(DAG$parents_previous_times))
  {
    var_idx$previous_times = lapply(DAG$parents_previous_times, function(idx)position_in_lower_tri_cross_vec(var_tag[idx], n_var = n_var, diag = T))
    var_idx$cross = mapply(
      function(x, y)
      {
        tatato = expand.grid(x, y)
        var_idx_12 = position_in_lower_tri(tatato[,1], tatato[,2], n_var)
      }
      ,
      lapply(mapply(c, DAG$parents_same_time, DAG$children, SIMPLIFY = F), function(x)var_tag[x]), 
      lapply(DAG$parents_previous_times                                 , function(x)var_tag[x]),# previous time first
      SIMPLIFY = F
    )
  }
  var_idx
}

##### lower_tri_idx_DAG = get_lower_tri_idx_DAG(DAG)
##### var_idx = get_var_idx(DAG, var_tag = var_tag)
##### # Pre-computing range and multiplier
##### multiplier = get_multiplier(
#####   a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
#####   r = r, A_vec = A_vec, nu_vec = nu_vec, a2_vec = a2_vec, u = u
##### )
##### effective_range = get_effective_range(
#####   a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
#####   r = r, A_vec = A_vec, a2_vec = a2_vec, u = u
##### )

get_linv_coeffs = function(
    locs_, var_tag, 
    DAG, 
    lower_tri_idx_DAG, 
    var_idx,
    multiplier, effective_range, 
    rho_vec, nu_vec
)
{
  res = list()
  res$coeffs = list()
  time_depth = nrow(multiplier)
  n_var = max(var_tag)
  nu_ = expand_nu(nu_vec)
  rho_vec_with_ones = put_ones_in_rho_vec(rho_vec)
  # looping over each block of the DAG
  #t1 = Sys.time()
  res = parallel::mclapply(
    X = seq(length(DAG[[1]])), mc.cores = parallel::detectCores(), 
    FUN = function(i)
    {
      print(i)
      covmat_coeffs_same_time = GMA_compressed(
        locs_ = locs_[c(DAG$parents_same_time[[i]], DAG$children[[i]]),, drop = F],
        lower_tri_idx = lower_tri_idx_DAG$same_time[[i]],
        var_idx = var_idx$current_time[[i]],
        multiplier = multiplier[1,,drop=FALSE], 
        effective_range = effective_range,
        nu_ = nu_,
        rho_vec_with_ones = rho_vec_with_ones
      )
      covmat_coeffs_previous = NULL
      covmat_coeffs_cross =  NULL
      if(nrow(multiplier)>1)
      {
        covmat_coeffs_previous = GMA_compressed(
          locs_ = locs_[DAG$parents_previous_times[[i]],, drop = F],
          lower_tri_idx = lower_tri_idx_DAG$previous_times[[i]],
          var_idx = var_idx$previous_times[[i]],
          multiplier = multiplier[-nrow(multiplier),,drop=FALSE], 
          effective_range = effective_range,
          nu_ = nu_,
          rho_vec_with_ones = rho_vec_with_ones
        )
        covmat_coeffs_cross =  
          GMA_rectangular(
            locs_1 = locs_[c(DAG$parents_same_time[[i]], DAG$children[[i]]),, drop = F], 
            locs_2 = locs_[DAG$parents_previous_times[[i]],, drop = F], 
            var_idx = var_idx$cross[[i]],
            multiplier = multiplier[-1,,drop=FALSE], 
            effective_range = effective_range[-1,,drop=FALSE],
            nu_ = expand_nu(nu_vec),
            n_var =  n_var, 
            rho_vec_with_ones = rho_vec_with_ones
          )
      }
      # my_image( expand_full_covmat(
      #   covmat_previous_periods = expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_previous$covmat,  block_lower_tri_idx = lower_tri_idx_DAG$previous_times[[i]], n_loc = length(DAG$parents_previous_times[[i]])),
      #   covmat_current_period   = expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_same_time$covmat, block_lower_tri_idx = lower_tri_idx_DAG$same_time     [[i]], n_loc = length(DAG$children[[i]])+length(DAG$parents_same_time[[i]])),
      #   side_blocks_rectangles  = covmat_coeffs_cross$covmat
      # ))
      covmat_chol = chol(
        expand_full_covmat(
          covmat_previous_periods = expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_previous$covmat,  block_lower_tri_idx = lower_tri_idx_DAG$previous_times[[i]], n_loc = length(DAG$parents_previous_times[[i]])),
          covmat_current_period   = expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_same_time$covmat, block_lower_tri_idx = lower_tri_idx_DAG$same_time     [[i]], n_loc = length(DAG$children[[i]])+length(DAG$parents_same_time[[i]])),
          side_blocks_rectangles  = covmat_coeffs_cross$covmat
        )
      )
      M = matrix(0, nrow(covmat_chol), length(DAG$children[[i]]))
      M[cbind(seq(nrow(covmat_chol)- length(DAG$children[[i]])+1, nrow(covmat_chol)), seq(length(DAG$children[[i]])))]=1
      return(backsolve(
        r =  covmat_chol, 
        M
      ))
      #my_image(covmat_chol)
      #### full_covmat = expand_full_covmat(
      ####   covmat_previous_periods = expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_previous$covmat,  block_lower_tri_idx = lower_tri_idx_DAG$previous_times[[i]], n_loc = length(DAG$parents_previous_times[[i]])),
      ####   covmat_current_period   = expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_same_time$covmat, block_lower_tri_idx = lower_tri_idx_DAG$same_time     [[i]], n_loc = length(DAG$children[[i]])+length(DAG$parents_same_time[[i]])),
      ####   side_blocks_rectangles  = 0*covmat_coeffs_cross$covmat
      #### )
      #### 1/sqrt(full_covmat[nrow(full_covmat), nrow(full_covmat)] - full_covmat[nrow(full_covmat), -nrow(full_covmat)]%*%solve(full_covmat[-nrow(full_covmat), -nrow(full_covmat)]) %*%full_covmat[nrow(full_covmat), -nrow(full_covmat)])
      #### (solve(full_covmat[-nrow(full_covmat), -nrow(full_covmat)]) %*%full_covmat[nrow(full_covmat), -nrow(full_covmat)])/c(sqrt(full_covmat[nrow(full_covmat), nrow(full_covmat)] - full_covmat[nrow(full_covmat), -nrow(full_covmat)]%*%solve(full_covmat[-nrow(full_covmat), -nrow(full_covmat)]) %*%full_covmat[nrow(full_covmat), -nrow(full_covmat)]))
      
      #plot(res$coeffs[[i]], solve(covmat_chol)[1,])
      # each row of solve covmat_chol factorizes the joint distribuion of children and parents
      # top rows = involves children + parents 
    }
  )
  #Sys.time()-t1
  res
}

##### Linv_coeffs = get_linv_coeffs(
#####   DAG = DAG, locs_ = locs_, 
#####   lower_tri_idx_DAG = lower_tri_idx_DAG, 
#####   var_idx = var_idx, 
#####   var_tag = var_tag, 
#####   multiplier = multiplier, 
#####   effective_range = effective_range, 
#####   rho_vec = rho_vec)
##### 
##### coeffs = Linv_coeffs

get_vecchia_blocks = function(DAG, coeffs, time_depth)
{
  same_time_indices = mapply(c, DAG$parents_same_time, DAG$children, SIMPLIFY = F)
  j = unlist(mapply(function(x, y)outer(x, rep(1, length(y))), same_time_indices, DAG$children))
  i = unlist(mapply(function(x, y)outer(rep(1, length(x)), y), same_time_indices, DAG$children))
  x = unlist(mapply(function(x, y)x[seq(nrow(x)-length(y)+1, nrow(x)),], coeffs, same_time_indices))# selecting bottom of  coeff matrices
  tatato = unlist(mapply(function(x, y)x[-seq(nrow(x)-length(y)+1, nrow(x)),], coeffs, same_time_indices))# selecting bottom of  coeff matrices
  triangular_on_diag = Matrix::sparseMatrix(
    i = i[i>=j],
    j = j[i>=j],
    x = x[i>=j],
    triangular = T
  )
  rectangular_below_diag = NULL
  if(time_depth > 1)
  {
    rectangular_below_diag = 
      lapply(seq(time_depth-1), function(idx)
        Matrix::sparseMatrix(
          j = unlist(mapply(function(x, y)outer(x, rep(1, length(y))), DAG$parents_previous_times, DAG$children)), 
          i = unlist(mapply(function(x, y)outer(rep(1, length(x)), y), DAG$parents_previous_times, DAG$children)),
          x = unlist(mapply(function(x, y)x[seq((idx-1)*length(y)+1,(idx)*length(y)),], coeffs, DAG$parents_previous_times, SIMPLIFY = F), recursive = T) 
        )
      )
    names(rectangular_below_diag)=sapply(seq(time_depth-1), function(i)paste("lag=", time_depth-i, sep=""))
  }
  list(
    triangular_on_diag = triangular_on_diag,
    rectangular_below_diag = rectangular_below_diag
  )
  
}  

# multiplication of a vector by the compressed repetitive sparse prior Chol
vecchia_blocks_mult = function(x, vecchia_blocks, time_depth)
{
  # space only case
  if(time_depth == 1)return(vecchia_blocks$triangular_on_diag %*% x)
  # putting vector in matrix
  x_ = matrix(x, nrow(vecchia_blocks$triangular_on_diag))
  # result matrix
  res = matrix(0, nrow(x_), ncol(x_))
  # first columns
  res[,1:(time_depth-1)]=x_[,1:(time_depth-1)]
  # multiplying current time blocks
  res[,time_depth:ncol(res)] = as.matrix(vecchia_blocks$triangular_on_diag %*% x_[,time_depth:ncol(res)])
  # multiplying previous time blocks
  for(i in seq(time_depth-1)) res[,time_depth:ncol(res)] = res[,time_depth:ncol(res)] + as.matrix(vecchia_blocks$rectangular_below_diag[[i]] %*% x_[,seq(0, ncol(res)-time_depth)+i])
  c(res)
}

# multiplication of a vector by the compressed repetitive sparse prior Chol
vecchia_blocks_t_mult = function(x, vecchia_blocks, time_depth)
{
  # putting vector in matrix
  x_ = matrix(x, nrow(vecchia_blocks$triangular_on_diag))
  # result matrix
  res = matrix(0, nrow(x_), ncol(x_))
  # first columns
  res[,1:(time_depth-1)]=x_[,1:(time_depth-1)]
  # multiplying current time blocks
  res[,time_depth:ncol(res)] = as.matrix(t(x_[,time_depth:ncol(res)]) %*% vecchia_blocks$triangular_on_diag)
  # multiplying previous time blocks
  for(i in seq(time_depth-1)) res[,time_depth:ncol(res)] = res[,time_depth:ncol(res)] + t(as.matrix(t(x_[,seq(0, ncol(res)-time_depth)+i]) %*% vecchia_blocks$rectangular_below_diag[[i]]))
  c(res)
}

# solving of a linear by the compressed repetitive sparse prior Chol
vecchia_blocks_solve = function(x, vecchia_blocks, time_depth)
{
  # space only case
  if(time_depth==1)return(as.vector(Matrix::solve(vecchia_blocks$triangular_on_diag, x)))
  # stacking x in  matrix
  x_ = matrix(x, nrow(vecchia_blocks$triangular_on_diag))
  # result
  res = matrix(0, nrow(x_), ncol(x_))
  # initializing result
  res[,1:(time_depth-1)]=x_[,1:(time_depth-1)]
  # solving recursively
  for(i in seq(time_depth, ncol(x_)))
  {
    res[,i] = x_[,i]
    # triangular block matrix solve
    for(j in seq(time_depth-1)) res[,i] = res[,i] - 
        as.matrix(vecchia_blocks$rectangular_below_diag[[j]] %*% res[,i-time_depth+j])
    res[,i] = as.vector(Matrix::solve(vecchia_blocks$triangular_on_diag, res[,i]))
  }
  c(res)
}



##### vecchia_blocks = get_vecchia_blocks(DAG, coeffs, time_depth)  
##### 
##### # checking t mult and mult
##### x1 = rnorm(100*n_loc)
##### x2 = rnorm(100*n_loc)
##### vecchia_blocks_t_mult(x1, vecchia_blocks, time_depth)*x2 - 
#####   vecchia_blocks_mult(x2, vecchia_blocks, time_depth)*x1
##### 
##### # checking solve and mult
##### x = rnorm(1000*n_loc)
##### x - vecchia_blocks_mult(vecchia_blocks_solve(x, vecchia_blocks, time_depth), vecchia_blocks, time_depth)
##### 
##### x = rnorm(100*n_loc)
##### w = vecchia_blocks_solve(x, vecchia_blocks, time_depth)




####  # Space only, no time
####  # flip image function to match matrix
####  my_image = function(m)image(t(m)[,nrow(m):1])
####  
####  # simulate data
####  set.seed(2)
####  
####  n_var = 3
####  time_depth = 1
####  
####  # spatial locations and variable index
####  locs_ = as.matrix(expand.grid(seq(100), seq(100)))/20#cbind(runif(n_loc), runif(n_loc))
####  locs_ = locs_ + rnorm(length(locs_), 0, .0001)
####  locs_ = locs_[sample(seq(nrow(locs_)), size = nrow(locs_), F),]
####  n_loc =nrow(locs_)
####  row.names(locs_) = seq(n_loc)
####  var_tag = ceiling(n_var*runif(n_loc));var_tag = match(var_tag, unique(var_tag))
####  
####  # covariance params
####  rho = GpGp::exponential_isotropic(c(1, 1, 0), .1*matrix(rnorm(2*n_var), n_var))
####  rho_vec = rho[lower.tri(rho, diag = F)];  remove(rho)
####  a2_vec = rep(100, n_var)#.000001*runif(n_var)
####  nu_vec = .5 + 2*runif(n_var)
####  alpha = .01 * runif(1)
####  a = runif(1) 
####  b = runif(1) 
####  cc = 1*runif(1) 
####  lambda = runif(1) 
####  delta = runif(1) 
####  r = 1 * runif(1) 
####  A_vec = rep(0, n_var)#runif(n_var)
####  u = seq(0, 0)
####  
####  # creating the DAG
####  NNarray_same_time = find_ordered_nn_multi(locs_ = locs_, 
####                                           var_tag = var_tag, 
####                                           m_whatever_closest = 10, m_same_var = 5, m_other_vars = 0, 
####                                           lonlat = F)
####  DAG = list(
####   "children" = lapply(seq(n_loc), function(x)x), 
####   "parents_same_time" = NNarray_same_time, 
####   "parents_previous_times" = 123456
####  )
####  
####  
####  lower_tri_idx_DAG = get_lower_tri_idx_DAG(DAG)
####  var_idx = get_var_idx(DAG, var_tag = var_tag)
####  # Pre-computing range and multiplier
####  multiplier = get_multiplier(
####    a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
####    r = r, A_vec = A_vec, nu_vec = nu_vec, a2_vec = a2_vec, u = u
####  )
####  effective_range = get_effective_range(
####    a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
####    r = r, A_vec = A_vec, a2_vec = a2_vec, u = u
####  )
####  
####  
####  Linv_coeffs = get_linv_coeffs(
####   DAG = DAG, locs_ = locs_, 
####   lower_tri_idx_DAG = lower_tri_idx_DAG, 
####   var_idx = var_idx, 
####   var_tag = var_tag, 
####   multiplier = multiplier, 
####   effective_range = effective_range, 
####   rho_vec = rho_vec)
#### ## 
#### coeffs = Linv_coeffs
#### vecchia_blocks = get_vecchia_blocks(DAG, coeffs, time_depth)  
####  
#### x = rnorm(10000)
#### x- vecchia_blocks_mult(vecchia_blocks_solve(x, vecchia_blocks = vecchia_blocks, time_depth = 1), vecchia_blocks = vecchia_blocks, time_depth = 1)
#### w = vecchia_blocks_solve(x, vecchia_blocks = vecchia_blocks, time_depth = 1)
#### Bidart::plot_pointillist_painting(locs_, w)
#### Bidart::plot_pointillist_painting(locs_[var_tag==1,], w[var_tag==1])
#### Bidart::plot_pointillist_painting(locs_[var_tag==2,], w[var_tag==2])
#### Bidart::plot_pointillist_painting(locs_[var_tag==3,], w[var_tag==3])


 



###  # Space only, no time in 1 dim
###  # flip image function to match matrix
###  my_image = function(m)image(t(m)[,nrow(m):1])
###  
###  # simulate data
###  set.seed(2)
###  
###  n_var = 8
###  time_depth = 1
###  
###  # spatial locations and variable index
###  locs_ = cbind(1, 10*runif(1000))
###  n_loc =nrow(locs_)
###  row.names(locs_) = seq(n_loc)
###  var_tag = ceiling(n_var*runif(n_loc));var_tag = match(var_tag, unique(var_tag))
###  
###  # covariance params
###  rho = GpGp::exponential_isotropic(c(1, 1, 0), .1*matrix(rnorm(2*n_var), n_var))
###  rho_vec = rho[lower.tri(rho, diag = F)];  remove(rho)
###  a2_vec = rep(1, n_var)#.000001*runif(n_var)
###  nu_vec = .5 + 2*runif(n_var)
###  alpha = .01 * runif(1)
###  a = runif(1) 
###  b = runif(1) 
###  cc = 1*runif(1) 
###  lambda = runif(1) 
###  delta = runif(1) 
###  r = 1 * runif(1) 
###  A_vec = rep(0, n_var)#runif(n_var)
###  u = seq(0, 0)
###  
###  # creating the DAG
###  NNarray_same_time = find_ordered_nn_multi(locs_ = locs_, 
###                                           var_tag = var_tag, 
###                                           m_whatever_closest = 10, m_same_var = 5, m_other_vars = 0, 
###                                           lonlat = F)
###  DAG = list(
###   "children" = lapply(seq(n_loc), function(x)x), 
###   "parents_same_time" = NNarray_same_time, 
###   "parents_previous_times" = 123456
###  )
###  
###  
###  lower_tri_idx_DAG = get_lower_tri_idx_DAG(DAG)
###  var_idx = get_var_idx(DAG, var_tag = var_tag)
###  # Pre-computing range and multiplier
###  multiplier = get_multiplier(
###    a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
###    r = r, A_vec = A_vec, nu_vec = nu_vec, a2_vec = a2_vec, u = u
###  )
###  effective_range = get_effective_range(
###    a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
###    r = r, A_vec = A_vec, a2_vec = a2_vec, u = u
###  )
###  
###  
###  Linv_coeffs = get_linv_coeffs(
###   DAG = DAG, locs_ = locs_, 
###   lower_tri_idx_DAG = lower_tri_idx_DAG, 
###   var_idx = var_idx, 
###   var_tag = var_tag, 
###   multiplier = multiplier, 
###   effective_range = effective_range, 
###   rho_vec = rho_vec)
### ## 
### coeffs = Linv_coeffs
### vecchia_blocks = get_vecchia_blocks(DAG, coeffs, time_depth)  
###  
### x = rnorm(n_loc)
### x- vecchia_blocks_mult(vecchia_blocks_solve(x, vecchia_blocks = vecchia_blocks, time_depth = 1), vecchia_blocks = vecchia_blocks, time_depth = 1)
### w = vecchia_blocks_solve(x, vecchia_blocks = vecchia_blocks, time_depth = 1)
### Bidart::plot_pointillist_painting(locs_, w)
### plot(locs_[var_tag==1,2], w[var_tag==1])
### plot(locs_[var_tag==2,2], w[var_tag==2])
### plot(locs_[var_tag==3,2], w[var_tag==3])
### plot(locs_[var_tag==4,2], w[var_tag==4])
### plot(locs_[var_tag==5,2], w[var_tag==5])
### plot(locs_[var_tag==6,2], w[var_tag==6])






 