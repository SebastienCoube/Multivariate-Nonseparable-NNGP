# Idea of field sampling: update all time periods of one cluster of spatial locations + the buffer at the beginning


# get precision blocks. Each block corresponds to a time lag, the first being lag = 0

get_precision_blocks = function(vecchia_blocks)
{
  res = list()
  res[[1]] = Matrix::crossprod(vecchia_blocks$triangular_on_diag) + Reduce("+", lapply(vecchia_blocks$rectangular_below_diag, Matrix::crossprod))
  k = length(vecchia_blocks$rectangular_below_diag)
  if(length(k)>0)
  {
    for(i in seq(k))
    {
      res[[i+1]] = Matrix::crossprod(vecchia_blocks$triangular_on_diag, vecchia_blocks$rectangular_below_diag[[i]])
      if(i!=k) res[[i+1]] = res[[i+1]] +  Reduce("+", mapply(Matrix::crossprod, vecchia_blocks$rectangular_below_diag[seq(k-i)], vecchia_blocks$rectangular_below_diag[seq(i+1, k)]))
    }
  }
  res
}


precision_blocks = get_precision_blocks(mcmc_nngp_list$chains[[1]]$stuff$vecchia)
noise_info = mcmc_nngp_list$chains[[1]]$stuff$noise_info
X_noise_list = mcmc_nngp_list$covariates$X_noise
diag_precision_block = precision_blocks[[1]]
useful_stuff = mcmc_nngp_list$useful_stuff

loc_subset_idx = seq(useful_stuff$n_loc/2)


get_a_posteriori_precision_chols = function(
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
  # time periods
  k = 1
  if(useful_stuff$n_time_periods>1) k = seq(useful_stuff$buffer_depth+1, useful_stuff$n_time_periods)
  a_postiori_pecision_chols = 
    lapply(k, 
           function(i)
           {
             posterior_precision = prior_precision
             #indices of non-na in y
             for(j in seq(length(loc_subset_idx)))
             {
               if(length(useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]]>0))
               {
                 tau_precision = matrix(0, n_field_per_site_subset[j], n_field_per_site_subset[j])
                 tau_precision[useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]],useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]]] = noise_info[[1]]
                 
                 cumsum_n_field_per_site_subset[j]+useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]]
               }
              }
             
           })
}

Vecchia_approx_DAG$field_position$location_idx 

# 


mcmc_nngp_list


clust = kmeans(locs[], centers = 1+floor(nrow(locs)/cluster_max_size))$clust





