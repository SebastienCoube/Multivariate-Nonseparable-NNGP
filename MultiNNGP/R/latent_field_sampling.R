



#posterior_precision_solve_chol : precision bits missing + check conformity with t solve


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

vecchia_blocks = mcmc_nngp_list$chains[[1]]$stuff$vecchia
precision_blocks = get_precision_blocks(vecchia_blocks)
noise_info = mcmc_nngp_list$chains[[1]]$stuff$noise_info
X_noise_list = mcmc_nngp_list$covariates$X_noise
diag_precision_block = precision_blocks[[1]]
useful_stuff = mcmc_nngp_list$useful_stuff

loc_subset_idx = seq(useful_stuff$n_loc/2)


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



get_precision_cholesky = function(
    precision_blocks, 
    noise_info, 
    X_noise_list, 
    useful_stuff, 
    loc_subset_idx, 
    time_begin, time_end, 
    Vecchia_approx_DAG
)
{
  a_posteriori_diag_precision = get_a_posteriori_diag_precision(
    diag_precision_block = precision_blocks[[1]],  
    noise_info = noise_info, X_noise_list = X_noise_list, 
    useful_stuff = useful_stuff, loc_subset_idx = loc_subset_idx, 
    Vecchia_approx_DAG = Vecchia_approx_DAG
  )
  selected_loc_var_pairs = unlist(Vecchia_approx_DAG$field_position$loc_match[loc_subset_idx]) 
  precision_blocks = lapply(precision_blocks, function(x)x[selected_loc_var_pairs,selected_loc_var_pairs])
  below_diag_blocks = precision_blocks; below_diag_blocks[[1]]=NULL
  # using block cholesky formula
  res = list()
  
  for(i in seq(length(a_posteriori_diag_precision)))
  {
    print("")
    print(paste("i=", i))
    print("")
    res[[i]] = list()
    # block diag
    if(i==1)res[[i]]$D = Matrix::expand(Matrix::Cholesky(a_posteriori_diag_precision[[1]]))
    if(i>1)
    {
      # coefficients above the block diag
      res[[i]]$L = list()
      cholesky_band_depth = min(i, useful_stuff$time_depth)-1
      below_diag_blocks_temp = below_diag_blocks[seq(cholesky_band_depth)]
      for(lag in seq(cholesky_band_depth, 1))
      {
        print(paste("lag=", lag))
        # recursive block solving starting by upper left corner
        # solving by upper left corner of previous cholesky transpose
        res[[i]]$L[[paste("lag=", lag, sep="")]] = 
          Matrix::t(res[[i-lag]]$D$P) %*%
          Matrix::solve(
            Matrix::t(res[[i-lag]]$D$L), 
            below_diag_blocks_temp[[lag]])
        # propagating solved block to blocks needing to be solved
        if(lag>1)
        {
          for(smaller_lag in seq(lag-1, 1))
          {
            below_diag_blocks_temp[[smaller_lag]] = 
              below_diag_blocks_temp[[smaller_lag]] -
              Matrix::crossprod(
                res[[i-smaller_lag]]$L[[lag-smaller_lag]], 
                res[[i]]$L[[paste("lag=", lag, sep="")]])
          }
        }
      }
      # solving last block, lower right
      res[[i]]$D = Matrix::expand(Matrix::Cholesky(
        a_posteriori_diag_precision[[i]] - Reduce("+", lapply(res[[i]]$L, Matrix::crossprod))
      ))
    }
  }
  
}







# part of th precision matrix off the diagonal, between the selected field and the observations corresponding to them
# In virtue of the data model, the rest of the observations is independent on the considered latent variables CONDITIONALLY ON the rest of the latent variables
get_field_obs_precision = function(
    noise_info, 
    X_noise_list, 
    useful_stuff, 
    loc_subset_idx, 
    Vecchia_approx_DAG,
    y
)
{
  #loc-var pairs of the DAG selected by loc only partition
  selected_loc_var_pairs = unlist(Vecchia_approx_DAG$field_position$loc_match[loc_subset_idx]) 
  # number of field variables simulated per site
  n_field_per_site_subset = unlist(useful_stuff$n_field_per_site[loc_subset_idx])
  cumsum_n_field_per_site_subset = cumsum(c(0, n_field_per_site_subset))
  # indices of additional precision due to observations
  additional_precision_indices = do.call(rbind, lapply(seq(length(loc_subset_idx)), function(j)
    {
    if(n_field_per_site_subset[j]==0)return(matrix(0, 0, 2))
     
      cbind(
        cumsum_n_field_per_site_subset[j] + c(outer(seq(   n_field_per_site_subset[j]), rep(1, dim(y)[2]))), 
        dim(y)[2]*(j-1) + c(outer(rep(1, n_field_per_site_subset[j]), seq(   dim(y)[2]))) 
      )
    }
    )
    )
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
            print(j)
            if(n_field_per_site_subset[j]==0)return(NULL)
            # creating matrix of 0 whose size corresponds to the number of latent field occurences at the selected location
            # filling the matrix with noise precision, at the indices corresponding to observed variables
            # eg: variables 1, 2, 4 are observed at site s. So tau_precision will have size 3*3. 
            # But at time t, only variables 1 and 4 are observed. So tau_precision will look like: 
            #    ?  0  ?
            #    0  0  0
            #    ?  0  ?
            tau_precision = matrix(0, n_field_per_site_subset[j], dim(y)[2])
            ##if(length(useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]])>1)
            #{
              #tau_precision[,useful_stuff$y_at_least_one_obs[[loc_subset_idx[j]]], drop = F][useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]],useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]]] = c(noise_info[[loc_subset_idx[[j]]]][[i]]$tau_precision)
              tau_precision[useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]] , seq(ncol(tau_precision))[useful_stuff$y_at_least_one_obs[[loc_subset_idx[j]]]] [useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]]] ] = c(noise_info[[loc_subset_idx[[j]]]][[i]]$tau_precision)
            #}
            
            c(tau_precision)
          }))
        Matrix::sparseMatrix(
          i = c(additional_precision_i),
          j = c(additional_precision_j), 
          x = c(additional_precision_x)
        ) 
      })
}




#get_a_posteriori_diag_precision_chols = function(
    #    diag_precision_block, 
#    noise_info, 
#    X_noise_list, 
#    useful_stuff, 
#    loc_subset_idx, 
#    Vecchia_approx_DAG
#)
#{
#  #loc-var pairs of the DAG selected by loc only partition
#  selected_loc_var_pairs = unlist(Vecchia_approx_DAG$field_position$loc_match[loc_subset_idx]) 
#  #subset of the all-loc-var-pairs precision
#  prior_precision = diag_precision_block[selected_loc_var_pairs,selected_loc_var_pairs]
#  # number of field variables simulated per site
#  n_field_per_site_subset = unlist(useful_stuff$n_field_per_site[loc_subset_idx])
#  cumsum_n_field_per_site_subset = cumsum(c(0, n_field_per_site_subset))
#  # indices of prior precision
#  prior_ij = sparseHessianFD::Matrix.to.Coord(prior_precision)
#  prior_i = prior_ij$rows[prior_ij$rows>=prior_ij$cols]
#  prior_j = prior_ij$cols[prior_ij$rows>=prior_ij$cols]
#  prior_x = prior_precision[cbind(prior_i, prior_j)]
#  # indices of additional precision due to observations
#  additional_precision_indices = do.call(
#    rbind, 
#    lapply(seq(length(loc_subset_idx)), function(j)
#      cumsum_n_field_per_site_subset[j] + 
#        get_lower_tri_idx(n_field_per_site_subset[j], diag = T)))
#  additional_precision_i = additional_precision_indices[,1]
#  additional_precision_j = additional_precision_indices[,2]
#  # time periods
#  k = 1
#  if(useful_stuff$n_time_periods>1) k = seq(useful_stuff$buffer_depth+1, useful_stuff$n_time_periods)
#  a_postiori_pecision_chols = 
#    parallel::mclapply(
#      mc.cores = parallel::detectCores()-1,
#      X = k, 
#      function(i)
#      {   
#        additional_precision_x = 
#          unlist(lapply(seq(length(loc_subset_idx)), function(j)
#          {
#            tau_precision = matrix(0, n_field_per_site_subset[j], n_field_per_site_subset[j])
#            if(length(useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]])>0)
#            {
#              # creating matrix of 0 whose size corresponds to the number of latent field occurences at the selected location
#              # filling the matrix with noise precision, at the indices corresponding to observed variables
#              # eg: variables 1, 2, 4 are observed at site s. So tau_precision will have size 3*3. 
#              # But at time t, only variables 1 and 4 are observed. So tau_precision will look like: 
#              #    ?  0  ?
#              #    0  0  0
#              #    ?  0  ?
#              tau_precision[useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]],useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]]] = noise_info[[loc_subset_idx[[j]]]][[i]]$tau_precision
#            }
#            tau_precision[lower.tri(tau_precision, diag = T)]
#          }))
#        Matrix::Cholesky(
#          Matrix::sparseMatrix(
#            symmetric = T,
#            i = c(prior_i, additional_precision_i),
#            j = c(prior_j, additional_precision_j), 
#            x = c(prior_x, additional_precision_x)
#          ) 
#        )
#      })
#}

# a_posteriori_diag_precision_chols = get_a_posteriori_diag_precision_chols(
#     diag_precision_block, 
#     noise_info, 
#     X_noise_list, 
#     useful_stuff, 
#     loc_subset_idx, 
#     Vecchia_approx_DAG
# )

## a_posteriori_diag_precision = get_a_posteriori_diag_precision(
##     diag_precision_block, 
##     noise_info, 
##     X_noise_list, 
##     useful_stuff, 
##     loc_subset_idx, 
##     Vecchia_approx_DAG
## )

# x = rnorm(mcmc_nngp_list$useful_stuff$buffer_depth*nrow(vecchia_blocks[[1]]) + nrow(a_posteriori_diag_precision_chols[[1]])*length(a_posteriori_diag_precision_chols))
# #1 time
# #x = rnorm(nrow(a_posteriori_diag_precision_chols[[1]])*length(a_posteriori_diag_precision_chols))
# below_diag_precision_blocks  =precision_blocks[-1]



