# Idea of field sampling: update all time periods of one cluster of spatial locations + the buffer at the beginning


# get precision blocks. Each block corresponds to a time lag, the first being lag = 0

get_precision_blocks = function(vecchia_blocks)
{
  if(length(vecchia_blocks$rectangular_below_diag)==0) return(list(Matrix::crossprod(vecchia_blocks$triangular_on_diag)))
  res = list()
  res[[1]] = Matrix::crossprod(vecchia_blocks$triangular_on_diag) + Reduce("+", lapply(vecchia_blocks$rectangular_below_diag, Matrix::crossprod))
  k = length(vecchia_blocks$rectangular_below_diag)
  for(i in seq(k))
  {
    res[[i+1]] = Matrix::crossprod(vecchia_blocks$triangular_on_diag, vecchia_blocks$rectangular_below_diag[[i]])
    if(i!=k) res[[i+1]] = res[[i+1]] +  Reduce("+", mapply(Matrix::crossprod, vecchia_blocks$rectangular_below_diag[seq(k-i)], vecchia_blocks$rectangular_below_diag[seq(i+1, k)]))
  }
  res
}

### vecchia_blocks = mcmc_nngp_list$chains[[1]]$stuff$vecchia
### precision_blocks = get_precision_blocks(vecchia_blocks)
### noise_info = mcmc_nngp_list$chains[[1]]$stuff$noise_info
### X_noise_list = mcmc_nngp_list$covariates$X_noise
### diag_precision_block = precision_blocks[[1]]
### useful_stuff = mcmc_nngp_list$useful_stuff
### 
### loc_subset_idx = seq(useful_stuff$n_loc/2)


get_a_posteriori_diag_precision = function(
    diag_precision_block, 
    noise_info, 
    X_noise_list, 
    useful_stuff, 
    loc_subset_idx, 
    Vecchia_approx_DAG
)
{
  #loc-var pairs of the DAG selected by loc only partition
  selected_loc_var_pairs = unlist(Vecchia_approx_DAG$field_position$loc_match[loc_subset_idx]) 
  #subset of the all-loc-var-pairs precision
  prior_precision = diag_precision_block[selected_loc_var_pairs,selected_loc_var_pairs]
  # number of field variables simulated per site
  n_field_per_site_subset = unlist(useful_stuff$n_field_per_site[loc_subset_idx])
  cumsum_n_field_per_site_subset = cumsum(c(0, n_field_per_site_subset))
  # indices of prior precision
  prior_ij = sparseHessianFD::Matrix.to.Coord(prior_precision)
  prior_i = prior_ij$rows[prior_ij$rows>=prior_ij$cols]
  prior_j = prior_ij$cols[prior_ij$rows>=prior_ij$cols]
  prior_x = prior_precision[cbind(prior_i, prior_j)]
  # indices of additional precision due to observations
  additional_precision_indices = do.call(rbind, lapply(seq(length(loc_subset_idx)), function(j)cumsum_n_field_per_site_subset[j] + get_lower_tri_idx(n_field_per_site_subset[j], diag = T)))
  additional_precision_i = additional_precision_indices[,1]
  additional_precision_j = additional_precision_indices[,2]
  # time periods
  k = 1
  if(useful_stuff$n_time_periods>1) k = seq(useful_stuff$buffer_depth+1, useful_stuff$n_time_periods)
  a_postiori_pecision_chols = 
    parallel::mclapply(
      mc.cores = parallel::detectCores()-1,
      X = k, 
      function(i)
      {   
        additional_precision_x = 
          unlist(lapply(seq(length(loc_subset_idx)), function(j)
          {
            tau_precision = matrix(0, n_field_per_site_subset[j], n_field_per_site_subset[j])
            if(length(useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]])>0)
            {
              # creating matrix of 0 whose size corresponds to the number of latent field occurences at the selected location
              # filling the matrix with noise precision, at the indices corresponding to observed variables
              # eg: variables 1, 2, 4 are observed at site s. So tau_precision will have size 3*3. 
              # But at time t, only variables 1 and 4 are observed. So tau_precision will look like: 
              #    ?  0  ?
              #    0  0  0
              #    ?  0  ?
              tau_precision[useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]],useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]]] = noise_info[[loc_subset_idx[[j]]]][[i]]$tau_precision
            }
            tau_precision[lower.tri(tau_precision, diag = T)]
          }))
        Matrix::sparseMatrix(
          symmetric = T,
          i = c(prior_i, additional_precision_i),
          j = c(prior_j, additional_precision_j), 
          x = c(prior_x, additional_precision_x)
        ) 
      })
}

get_a_posteriori_diag_precision_chols = function(
    diag_precision_block, 
    noise_info, 
    X_noise_list, 
    useful_stuff, 
    loc_subset_idx, 
    Vecchia_approx_DAG
    )
{
  #loc-var pairs of the DAG selected by loc only partition
  selected_loc_var_pairs = unlist(Vecchia_approx_DAG$field_position$loc_match[loc_subset_idx]) 
  #subset of the all-loc-var-pairs precision
  prior_precision = diag_precision_block[selected_loc_var_pairs,selected_loc_var_pairs]
  # number of field variables simulated per site
  n_field_per_site_subset = unlist(useful_stuff$n_field_per_site[loc_subset_idx])
  cumsum_n_field_per_site_subset = cumsum(c(0, n_field_per_site_subset))
  # indices of prior precision
  prior_ij = sparseHessianFD::Matrix.to.Coord(prior_precision)
  prior_i = prior_ij$rows[prior_ij$rows>=prior_ij$cols]
  prior_j = prior_ij$cols[prior_ij$rows>=prior_ij$cols]
  prior_x = prior_precision[cbind(prior_i, prior_j)]
  # indices of additional precision due to observations
  additional_precision_indices = do.call(rbind, lapply(seq(length(loc_subset_idx)), function(j)cumsum_n_field_per_site_subset[j] + get_lower_tri_idx(n_field_per_site_subset[j], diag = T)))
  additional_precision_i = additional_precision_indices[,1]
  additional_precision_j = additional_precision_indices[,2]
  # time periods
  k = 1
  if(useful_stuff$n_time_periods>1) k = seq(useful_stuff$buffer_depth+1, useful_stuff$n_time_periods)
  a_postiori_pecision_chols = 
    parallel::mclapply(
      mc.cores = parallel::detectCores()-1,
      X = k, 
      function(i)
      {   
        additional_precision_x = 
          unlist(lapply(seq(length(loc_subset_idx)), function(j)
          {
            tau_precision = matrix(0, n_field_per_site_subset[j], n_field_per_site_subset[j])
            if(length(useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]])>0)
            {
              # creating matrix of 0 whose size corresponds to the number of latent field occurences at the selected location
              # filling the matrix with noise precision, at the indices corresponding to observed variables
              # eg: variables 1, 2, 4 are observed at site s. So tau_precision will have size 3*3. 
              # But at time t, only variables 1 and 4 are observed. So tau_precision will look like: 
              #    ?  0  ?
              #    0  0  0
              #    ?  0  ?
              tau_precision[useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]],useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]]] = noise_info[[loc_subset_idx[[j]]]][[i]]$tau_precision
            }
            tau_precision[lower.tri(tau_precision, diag = T)]
          }))
        Matrix::Cholesky(
          Matrix::sparseMatrix(
            symmetric = T,
            i = c(prior_i, additional_precision_i),
            j = c(prior_j, additional_precision_j), 
            x = c(prior_x, additional_precision_x)
          ) 
        )
      })
}


#### test for block diagonal addition
###tatato = #parallel::mc
###lapply(
###  #mc.cores = parallel::detectCores()-1,
###  X = k, 
###  function(i)
###  {
###              prior_precision +
###                Matrix::bdiag(
###    lapply(seq(length(loc_subset_idx)), function(j)
###    {
###      tau_precision = matrix(0, n_field_per_site_subset[j], n_field_per_site_subset[j])
###      if(length(useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]])>0)
###      {
###        # creating matrix of 0 whose size corresponds to the number of latent field occurences at the selected location
###        # filling the matrix with noise precision, at the indices corresponding to observed variables
###        # eg: variables 1, 2, 4 are observed at site s. So tau_precision will have size 3*3. 
###        # But at time t, only variables 1 and 4 are observed. So tau_precision will look like: 
###        #    ?  0  ?
###        #    0  0  0
###        #    ?  0  ?
###        tau_precision[useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]],useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]]] = noise_info[[loc_subset_idx[[j]]]][[i]]$tau_precision
###        tau_precision
###      }
###    })
###    )
###  })

#### a_posteriori_diag_precision_chols = get_a_posteriori_diag_precision_chols(
####     diag_precision_block, 
####     noise_info, 
####     X_noise_list, 
####     useful_stuff, 
####     loc_subset_idx, 
####     Vecchia_approx_DAG
#### )

#### a_posteriori_diag_precision = get_a_posteriori_diag_precision(
####     diag_precision_block, 
####     noise_info, 
####     X_noise_list, 
####     useful_stuff, 
####     loc_subset_idx, 
####     Vecchia_approx_DAG
#### )

### x = rnorm(mcmc_nngp_list$useful_stuff$buffer_depth*nrow(vecchia_blocks[[1]]) + nrow(a_posteriori_diag_precision_chols[[1]])*length(a_posteriori_diag_precision_chols))
### #1 time
### #x = rnorm(nrow(a_posteriori_diag_precision_chols[[1]])*length(a_posteriori_diag_precision_chols))
### below_diag_precision_blocks  =precision_blocks[-1]


posterior_precision_solve_chol = function(
    a_posteriori_diag_precision_chols, # a posteriori blocks on the diagonal
    below_diag_precision_blocks,  # a priori blocks below the diagonal
    vecchia_blocks, # vecchia prior to solve the buffer
    loc_subset_idx, # subset of locations
    x # to be solved
){
  # space only case
  if(length(below_diag_precision_blocks)==0)return(as.vector(
    Matrix::t(as(a_posteriori_diag_precision_chols[[1]], "pMatrix")) %*% 
      Matrix::solve(
        Matrix::t(as(a_posteriori_diag_precision_chols[[1]], "Matrix")), 
        x
      )))
  # loc var pairs corresponding to subset loc
  selected_loc_var_pairs = unlist(Vecchia_approx_DAG$field_position$loc_match[loc_subset_idx]) 
  below_diag_precision_blocks_subset = lapply(below_diag_precision_blocks, function(x)x[selected_loc_var_pairs,selected_loc_var_pairs])
  below_diag_precision_blocks_cross = lapply(below_diag_precision_blocks, function(x)x[,selected_loc_var_pairs])
  
  x_buffer   = x[ seq(length(x) - nrow(a_posteriori_diag_precision_chols[[1]])*length(a_posteriori_diag_precision_chols))]
  x_observed = x[-seq(length(x) - nrow(a_posteriori_diag_precision_chols[[1]])*length(a_posteriori_diag_precision_chols))]
  x_buffer   = matrix(x_buffer,   nrow(below_diag_precision_blocks[[1]]))
  x_observed = matrix(x_observed, nrow(a_posteriori_diag_precision_chols[[1]]))
  # begining the solving with the buffer
  x_buffer = vecchia_blocks_solve(x_buffer, vecchia_blocks = vecchia_blocks, time_depth = length(precision_blocks)-1)
  # using block solve formula, starting from top
  for(i in seq(ncol(x_observed))){
    # multiplying previously computed elements by blocks above the diag
    for(j in seq(length(below_diag_precision_blocks_subset)))
    # solving by lower tri blocks on the diag
    x_observed[,i] = as.vector(
      Matrix::t(as(a_posteriori_diag_precision_chols[[i]], "pMatrix")) %*%
      Matrix::solve(
      Matrix::t(as(a_posteriori_diag_precision_chols[[i]], "Matrix")), 
      x_observed[,i]
      ))
  }
  c(x_buffer, x_observed)
}

posterior_precision_solve_t_chol = function(
    a_posteriori_diag_precision_chols, # a posteriori blocks on the diagonal
    below_diag_precision_blocks,  # a priori blocks below the diagonal
    vecchia_blocks, # vecchia prior to solve the buffer
    loc_subset_idx, # subset of locations
    x # to be solved
){
  if(length(below_diag_precision_blocks)==0)return(as.vector(Matrix::solve(
    as(a_posteriori_diag_precision_chols[[1]], "Matrix"), 
    as(a_posteriori_diag_precision_chols[[1]], "pMatrix") %*% 
      x
  )))
  # loc var pairs corresponding to subset loc
  selected_loc_var_pairs = unlist(Vecchia_approx_DAG$field_position$loc_match[loc_subset_idx]) 
  below_diag_precision_blocks_subset = lapply(below_diag_precision_blocks, function(x)x[selected_loc_var_pairs,selected_loc_var_pairs])
  below_diag_precision_blocks_cross = lapply(below_diag_precision_blocks, function(x)x[,selected_loc_var_pairs])
  
  x_buffer   = x[ seq(length(x) - nrow(a_posteriori_diag_precision_chols[[1]])*length(a_posteriori_diag_precision_chols))]
  x_observed = x[-seq(length(x) - nrow(a_posteriori_diag_precision_chols[[1]])*length(a_posteriori_diag_precision_chols))]
  x_buffer   = matrix(x_buffer,   nrow(below_diag_precision_blocks[[1]]))
  x_observed = matrix(x_observed, nrow(a_posteriori_diag_precision_chols[[1]]))
  # using block solve formula, starting from bottom
  for(i in seq(ncol(x_observed))){
    # multiplying previously computed elements by blocks above the diag
    if(i>1) for(j in seq(min(length(below_diag_precision_blocks), i-1)))x_observed[,ncol(x_observed)-i+1] = as.vector(x_observed[,ncol(x_observed)-i+1] -below_diag_precision_blocks_subset[[j]] %*% x_observed[,ncol(x_observed)-i+1+j]) 
    # solving by lower tri blocks on the diag
    x_observed[,ncol(x_observed)-i+1] = as.vector(Matrix::solve(
      as(a_posteriori_diag_precision_chols[[ncol(x_observed)-i+1]], "Matrix"), 
      as(a_posteriori_diag_precision_chols[[ncol(x_observed)-i+1]], "pMatrix") %*% 
      x_observed[,ncol(x_observed)-i+1]
      ))
  }
 # finishing the solving with the buffer
 for(j in seq(min(length(below_diag_precision_blocks), ncol(x_observed))))
 {
   x_buffer[,seq(ncol(x_buffer)-length(below_diag_precision_blocks)+j, ncol(x_buffer))] = 
   as.matrix(x_buffer[,seq(ncol(x_buffer)-length(below_diag_precision_blocks)+j, ncol(x_buffer))] + 
     below_diag_precision_blocks_cross [[j]] %*% x_observed[,seq(length(below_diag_precision_blocks)-j+1)])
 }
  x_buffer = vecchia_blocks_solve(x_buffer, vecchia_blocks = vecchia_blocks, time_depth = length(precision_blocks)-1)
  c(x_buffer, x_observed)
}

## one time test
#a_posteriori_diag_precision[[1]]%*%
#posterior_precision_solve_chol(
#  a_posteriori_diag_precision_chols = a_posteriori_diag_precision_chols, 
#  below_diag_precision_blocks = below_diag_precision_blocks, 
#  vecchia_blocks = vecchia_blocks, 
#  loc_subset_idx = loc_subset_idx, 
#  x = posterior_precision_solve_t_chol(
#    a_posteriori_diag_precision_chols = a_posteriori_diag_precision_chols, 
#    below_diag_precision_blocks = below_diag_precision_blocks, 
#    vecchia_blocks = vecchia_blocks, 
#    loc_subset_idx = loc_subset_idx, 
#    x = x)) - x

## checking that solve_chol and solve_t_chol do the same thing
#x1 = rnorm(mcmc_nngp_list$useful_stuff$buffer_depth*nrow(vecchia_blocks[[1]]) + nrow(a_posteriori_diag_precision_chols[[1]])*length(a_posteriori_diag_precision_chols))
#x2 = rnorm(mcmc_nngp_list$useful_stuff$buffer_depth*nrow(vecchia_blocks[[1]]) + nrow(a_posteriori_diag_precision_chols[[1]])*length(a_posteriori_diag_precision_chols))
#posterior_precision_solve_chol(
#a_posteriori_diag_precision_chols = a_posteriori_diag_precision_chols, 
#below_diag_precision_blocks = below_diag_precision_blocks, 
#vecchia_blocks = vecchia_blocks, 
#loc_subset_idx = loc_subset_idx, 
#x = x2)*
#posterior_precision_solve_t_chol(
#a_posteriori_diag_precision_chols = a_posteriori_diag_precision_chols, 
#below_diag_precision_blocks = below_diag_precision_blocks, 
#vecchia_blocks = vecchia_blocks, 
#loc_subset_idx = loc_subset_idx, 
#x = x1)-
#posterior_precision_solve_chol(
#a_posteriori_diag_precision_chols = a_posteriori_diag_precision_chols, 
#below_diag_precision_blocks = below_diag_precision_blocks, 
#vecchia_blocks = vecchia_blocks, 
#loc_subset_idx = loc_subset_idx, 
#x = x1)*
#posterior_precision_solve_t_chol(
#a_posteriori_diag_precision_chols = a_posteriori_diag_precision_chols, 
#below_diag_precision_blocks = below_diag_precision_blocks, 
#vecchia_blocks = vecchia_blocks, 
#loc_subset_idx = loc_subset_idx, 
#x = x2)



# part of th precision matrix off the diagonal, between the selected field and the observations corresponding to them
# In virtue of the data model, the rest of the observations is independent on the considered latent variables CONDITIONALLY ON the rest of the latent variables
get_field_obs_precision = function(
    noise_info, 
    X_noise_list, 
    useful_stuff, 
    loc_subset_idx, 
    Vecchia_approx_DAG
)
{
  #loc-var pairs of the DAG selected by loc only partition
  selected_loc_var_pairs = unlist(Vecchia_approx_DAG$field_position$loc_match[loc_subset_idx]) 
  # number of field variables simulated per site
  n_field_per_site_subset = unlist(useful_stuff$n_field_per_site[loc_subset_idx])
  cumsum_n_field_per_site_subset = cumsum(c(0, n_field_per_site_subset))
  # indices of additional precision due to observations
  additional_precision_indices = do.call(rbind, lapply(seq(length(loc_subset_idx)), function(j)cumsum_n_field_per_site_subset[j] + get_lower_tri_idx(n_field_per_site_subset[j], diag = T)))
  additional_precision_i = additional_precision_indices[,1]
  additional_precision_j = additional_precision_indices[,2]
  # time periods
  k = 1
  if(useful_stuff$n_time_periods>1) k = seq(useful_stuff$buffer_depth+1, useful_stuff$n_time_periods)
  field_obs_precision = 
    parallel::mclapply(
      mc.cores = parallel::detectCores()-1,
      X = k, 
      function(i)
      {   
        additional_precision_x = 
          unlist(lapply(seq(length(loc_subset_idx)), function(j)
          {
            tau_precision = matrix(0, n_field_per_site_subset[j], n_field_per_site_subset[j])
            if(length(useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]])>0)
            {
              # creating matrix of 0 whose size corresponds to the number of latent field occurences at the selected location
              # filling the matrix with noise precision, at the indices corresponding to observed variables
              # eg: variables 1, 2, 4 are observed at site s. So tau_precision will have size 3*3. 
              # But at time t, only variables 1 and 4 are observed. So tau_precision will look like: 
              #    ?  0  ?
              #    0  0  0
              #    ?  0  ?
              tau_precision[useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]],useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]]] = noise_info[[loc_subset_idx[[j]]]][[i]]$tau_precision
            }
            tau_precision[lower.tri(tau_precision, diag = T)]
          }))
        Matrix::sparseMatrix(
          symmetric = T,
          i = c(additional_precision_i),
          j = c(additional_precision_j), 
          x = c(additional_precision_x)
        ) 
      })
}


