
# useful_stuff = mcmc_nngp_list$useful_stuff
# noise_info = mcmc_nngp_list$chains$chain_1$stuff$noise_info
# X_noise_list = mcmc_nngp_list$covariates$X_noise



# multiplies z by the block-diagonal precision between the latent field and the observations
# z is under the same format as useful_stuff$y_split  , that is a matrix of lists
field_obs_precision_mult = function(
    noise_info, 
    X_noise_list, 
    useful_stuff, 
    Vecchia_approx_DAG, 
    z
)
{
  cumsum_n_field_per_site = cumsum(c(0, useful_stuff$n_field_per_site))
  # number of field variables simulated per site
  res_ = do.call(
    cbind, 
    parallel::mclapply(
    mc.cores = 10, 
    X = seq(useful_stuff$n_time_periods), 
    function(i_time)
    {
      res = rep(0, useful_stuff$n_field)
      for(i_site in seq(useful_stuff$n_loc))
      {
       if(length(useful_stuff$position_in_y_at_least_one_obs[[i_time]][[i_site]])>0)
       {
         idx =  useful_stuff$position_in_y_at_least_one_obs[[i_time]][[i_site]]
         res[seq(cumsum_n_field_per_site[i_site]+1, cumsum_n_field_per_site[i_site+1])[idx]] = 
           noise_info[[i_site]][[i_time]]$tau_precision %*% na.omit(z[i_site, i_time][[1]])
       }
      }
      res
    }
  ))
    
  ###res =   matrix(0, useful_stuff$n_field, useful_stuff$n_time_periods)
  ###for(i_site in seq(useful_stuff$n_loc))
  ###{
  ###  for(i_time in seq(useful_stuff$n_time_periods))
  ###  {
  ###    if(length(useful_stuff$position_in_y_at_least_one_obs[[i_time]][[i_site]])>0)
  ###    {
  ###      idx =  useful_stuff$position_in_y_at_least_one_obs[[i_time]][[i_site]]
  ###      res[seq(cumsum_n_field_per_site[i_site]+1, cumsum_n_field_per_site[i_site+1])[idx], i_time] = 
  ###        noise_info[[i_site]][[i_time]]$tau_precision %*% na.omit(useful_stuff$y_split[i_site, i_time][[1]])
  ###    }
  ###  }
  ###}
  ###summary(
  ###c(res)-
  ###c(res_)
  ###)
  res_
}

z = mcmc_nngp_list$y
field_obs_precision_mult(
    noise_info, 
    X_noise_list, 
    useful_stuff, 
    Vecchia_approx_DAG, 
    z
)

# a priori whitened gradient, that is design matrix equal to prior covariance

latent_field_grad = function(
    mcmc_nngp_list, 
    field, split_latent_field,
    vecchia_blocks
) {
    # Q = prior precision
    # N = noise precision
    # density = -.5 wT Q w -.5 (w-obs)^T N (w-obs) 
    # gradient = -.5 Q w -.5 N (w-obs) 
    # whitened gradient = -.5 w -.5 Q-1 N (w-obs) 
    
    # design matrix = a
  # gradient
  (  
    -.5 * field 
    + vecchia_blocks_solve(
      vecchia_blocks_t_solve(
        latent_field_data_gradient(split_latent_field, chain, mcmc_nngp_list)
        , vecchia_blocks = vecchia_blocks)
      , vecchia_blocks = vecchia_blocks
    )
  )
}

#y_minus_field = function(
#    field, mcmc_nngp_list
#    )
#{
#  res = mcmc_nngp_list$y
#  res[]=0
#  for(i in seq(ncol(field)))
#  {
#    y[mcmc_nngp_list$Vecchia_approx_DAG$field_position$location_idx, mcmc_nngp_list$Vecchia_approx_DAG$field_position$var_idx,] =
#    y[mcmc_nngp_list$Vecchia_approx_DAG$field_position$location_idx, mcmc_nngp_list$Vecchia_approx_DAG$field_position$var_idx,] - 
#      field[,i]
#  }
#}
#y_minus_field(proposed_field, mcmc_nngp_list)


# splits latent field by locations (following the column)
# and time (following the row)
get_split_latent_field = function(mcmc_nngp_list, field)
{
  lapply(tapply( split(field, col(field)), mcmc_nngp_list$Vecchia_approx_DAG$field_position$location_idx, c), function(x)do.call(rbind, x))
}


#split_latent_field = get_split_latent_field(mcmc_nngp_list, mcmc_nngp_list$chains$chain_1$params$field)

# gets the gradient of the density of the Gaussian observations with respect to the latent field 
latent_field_data_gradient = function(split_latent_field, chain, mcmc_nngp_list)
{
  #do.call(cbind, 
  unlist(
    parallel::mclapply(mc.cores = 10, 
                       seq(mcmc_nngp_list$useful_stuff$n_loc), function(i_loc)
                       {
                         # field - obs
                         z = lapply( 
                           split(
                             split_latent_field[[i_loc]]-(mcmc_nngp_list$y[i_loc,mcmc_nngp_list$useful_stuff$y_at_least_one_obs[[i_loc]],]), 
                             col(split_latent_field[[i_loc]])
                           ), 
                           na.omit
                         )
                         do.call(rbind, lapply(
                           seq(mcmc_nngp_list$useful_stuff$n_time_periods), 
                           function(i_time)
                           {
                             res = rep(0, mcmc_nngp_list$useful_stuff$n_field_per_site[[i_loc]])
                             if(!is.null(chain$stuff$noise_info[[i_loc]][[i_time]]$tau_precision)) res[mcmc_nngp_list$useful_stuff$position_in_y_at_least_one_obs[[i_time]][[i_loc]]] = res[mcmc_nngp_list$useful_stuff$position_in_y_at_least_one_obs[[i_time]][[i_loc]]] -.5 * z[[i_time]]%*% chain$stuff$noise_info[[i_loc]][[i_time]]$tau_precision
                             res
                           }
                         ))
      # y - latent field
  }))
}

#split_latent_field = get_split_latent_field(mcmc_nngp_list, chain)
latent_field_data_density = function(split_latent_field, chain, mcmc_nngp_list)
{
  sum(unlist(parallel::mclapply(mc.cores = 10, seq(mcmc_nngp_list$useful_stuff$n_loc), function(i_loc)
  {
   sum(
      mapply(
      function(z, precision)
      {
        if(is.null(precision))return(0)
        -.5 * sum(z * (z%*%precision))
      },
      # y - latent field
      z = lapply( 
             split(
               split_latent_field[[i_loc]]-(mcmc_nngp_list$y[i_loc,mcmc_nngp_list$useful_stuff$y_at_least_one_obs[[i_loc]],]), 
               col(split_latent_field[[i_loc]])
             ), 
             na.omit
      ), 
      precision = lapply(chain$stuff$noise_info[[i_loc]], function(x)x$tau_precision)
    ))
  })))
}

