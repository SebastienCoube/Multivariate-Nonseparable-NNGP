
homeostasy  = function(log_var, accepted, target  = .25, eps = .01)
{
  if(accepted) log_var = log_var + eps * (1-target)
  if(!accepted) log_var = log_var - eps * target
  return(log_var)
}


list2env(
  get_blocks_n_bases_with_witnesses(
    mcmc_nngp_list = mcmc_nngp_list, 
    space_cluster_size_target = 120, 
    chain = chain, which_vars = seq(6)),
  envir =  environment()
)
list2env(
  get_time_split(
    time_depth = mcmc_nngp_list$useful_stuff$time_depth, 
    n_time_periods = mcmc_nngp_list$useful_stuff$n_time_periods, 
    time_target_size = 100),
  envir =  environment()
)
t1 = Sys.time()
tatato = latent_field_blocked_sampling(
  field = chain$params$field, epsilon = NULL, 
  mcmc_nngp_list = mcmc_nngp_list, 
  time_coloring = time_coloring, spatial_coloring = spatial_coloring, 
  time_split = time_split, basis_functions = basis_functions, 
  precision_blocks = chain$stuff$precision_blocks, vecchia_blocks = chain$stuff$vecchia_blocks, 
  transposed_vecchia_blocks = chain$stuff$transposed_vecchia_blocks, noise_precisions = chain$stuff$noise_precisions)
print(Sys.time()-t1)

plot(rep(stuff_for_plots$locs_no_na[,1], stuff_for_plots$n_var), stuff_for_plots$y_true[,,80], 
     col = rep(c("lightgray", "lightpink", "lightgreen", "lightblue", "lightcyan", "lightpink"), each = stuff_for_plots$n_loc), 
     cex = .3, pch = 15, 
     xlab = "spatial site", 
     ylab = "latent field and true field"
)
points(mcmc_nngp_list$locs[mcmc_nngp_list$Vecchia_approx_DAG$field_position$location_idx,1], 
       tatato$field[,mcmc_nngp_list$useful_stuff$buffer_depth + 80], 
       cex = .3, pch = 1, col = mcmc_nngp_list$Vecchia_approx_DAG$field_position$var_idx
)



chain$params$field = 0 * chain$params$field





do_100_updates = function(chain, mcmc_nngp_list, kernel_learning_rate, thinning_rate = 1)
{
  iter = 1
  records = lapply(chain$params, function(x)array(0, c(dim(x), 200)))
  for(iter in seq(200))
  {
    
    #########################
    # Covariance parameters #
    #########################
    var_idx = 1
    for(var_idx in seq(mcmc_nngp_list$useful_stuff$n_var_y))
    {
      print(paste("var_idx=", var_idx))
      
      # sampling proposal ####
      proposal = rnorm(2, sd = exp(.5*chain$kernels$var_wise[var_idx]))
      # proposing range ####
      proposed_log_range_vec = chain$params$log_range_vec
#      proposed_log_range_vec[var_idx] = 
#        proposed_log_range_vec[var_idx] +
#        (
#          mcmc_nngp_list$hierarchical_model$log_range_prior[var_idx, 2] - 
#            mcmc_nngp_list$hierarchical_model$log_range_prior[var_idx, 1]
#        )* proposal[1]
      # proposing smoothness ####
      proposed_smoothness_vec = chain$params$smoothness_vec
#      proposed_smoothness_vec[var_idx] = 
#        proposed_smoothness_vec[var_idx] +
#        (
#          mcmc_nngp_list$hierarchical_model$smoothness_prior[var_idx, 2] - 
#            mcmc_nngp_list$hierarchical_model$smoothness_prior[var_idx, 1]
#        )* proposal[2]
      
      # checking prior ####
      accepted = F
      if(
        #####
        (proposed_log_range_vec[var_idx]>mcmc_nngp_list$hierarchical_model$log_range_prior[var_idx, 1])&
        (proposed_log_range_vec[var_idx]<mcmc_nngp_list$hierarchical_model$log_range_prior[var_idx, 2])&
        (proposed_smoothness_vec[var_idx]>mcmc_nngp_list$hierarchical_model$smoothness_prior[var_idx, 1])&
        (proposed_smoothness_vec[var_idx]<mcmc_nngp_list$hierarchical_model$smoothness_prior[var_idx, 2])
        #####
      ){
        #print("prior accepted")
        # proposing Vecchia factors ####
        proposed_vecchia_blocks = 
          vecchia_block_approx( 
            Vecchia_approx_DAG = mcmc_nngp_list$Vecchia_approx_DAG, locs = mcmc_nngp_list$locs, 
            time_depth = mcmc_nngp_list$useful_stuff$time_depth, #does not depend on params
            rho_vec = chain$params$rho_vec, a =  chain$params$a_scal,
            b = chain$params$b_scal, cc = chain$params$cc, delta = chain$params$delta, 
            lambda = chain$params$lambda, r = chain$params$r, 
            A_vec = chain$params$A_vec, 
            nu_vec = proposed_smoothness_vec, 
            log_range_vec = proposed_log_range_vec
          )
        proposed_transposed_vecchia_blocks = vecchia_blocks_t(proposed_vecchia_blocks)
        proposed_precision_blocks = get_precision_blocks(proposed_vecchia_blocks)
        # spatial and temporal saucissoning ##### 
        list2env(
          get_blocks_n_bases_with_witnesses(
            mcmc_nngp_list = mcmc_nngp_list, 
            space_cluster_size_target = 3, 
            chain = chain, which_vars = var_idx),
          envir =  environment()
        )
#basis_functions = list(Matrix::sparseMatrix(
#  i = which(Vecchia_approx_DAG$field_position$var_idx==1), 
#  j = seq(sum(Vecchia_approx_DAG$field_position$var_idx==1)), 
#  dims = c(mcmc_nngp_list$useful_stuff$n_field, sum(Vecchia_approx_DAG$field_position$var_idx==1))
#  ))
## basis_functions = basis_functions [1]
## spatial_coloring=1
        list2env(
          get_time_split(
            time_depth = mcmc_nngp_list$useful_stuff$time_depth, 
            n_time_periods = mcmc_nngp_list$useful_stuff$n_time_periods, 
            time_target_size = 100),
          envir =  environment()
        )
#time_split[1,2] = 25
#time_split[1,1] = 5
        # proposal of new field and transition ratio #####
        field_proposal = 
          latent_field_blocked_sampling(
            field = chain$params$field, epsilon = NULL, 
            mcmc_nngp_list = mcmc_nngp_list, 
            time_coloring = time_coloring, spatial_coloring = spatial_coloring, 
            time_split = time_split, basis_functions = basis_functions, 
            precision_blocks = proposed_precision_blocks, vecchia_blocks = proposed_vecchia_blocks, 
            transposed_vecchia_blocks = proposed_transposed_vecchia_blocks, 
            noise_precisions = chain$stuff$noise_precisions)
        reverse_field_proposal_density = 
          latent_field_blocked_sampling(
            field = field_proposal$field, epsilon = field_proposal$epsilon, 
            mcmc_nngp_list = mcmc_nngp_list, 
            time_coloring = time_coloring, spatial_coloring = spatial_coloring, 
            time_split = time_split, basis_functions = basis_functions, 
            precision_blocks = chain$stuff$precision_blocks, vecchia_blocks = chain$stuff$vecchia_blocks, 
            transposed_vecchia_blocks = chain$stuff$transposed_vecchia_blocks, 
            noise_precisions = chain$stuff$noise_precisions)
        # checking the reverse field proposal goes back to initial field
        # summary(c(reverse_field_proposal_density$field - chain$params$field))
        start_density = covparms_and_latent_field_density(
          field = chain$params$field, 
          vecchia_blocks = chain$stuff$vecchia_blocks, 
          noise_precisions = chain$stuff$noise_precisions, 
          mcmc_nngp_list = mcmc_nngp_list
        )
        end_density = covparms_and_latent_field_density(
          field = field_proposal$field, 
          vecchia_blocks = proposed_vecchia_blocks, 
          noise_precisions = chain$stuff$noise_precisions, 
          mcmc_nngp_list = mcmc_nngp_list
        )
        
        plot(rep(stuff_for_plots$locs_no_na[,1], stuff_for_plots$n_var), stuff_for_plots$y_true[,,10], 
            col = rep(c("lightgray", "lightpink", "lightgreen", "lightblue", "lightcyan", "lightpink"), each = stuff_for_plots$n_loc), 
            cex = .3, pch = 15, 
            xlab = "spatial site", 
            ylab = "latent field and true field"
        )
        points(mcmc_nngp_list$locs[mcmc_nngp_list$Vecchia_approx_DAG$field_position$location_idx,1], 
              field_proposal$field[,mcmc_nngp_list$useful_stuff$buffer_depth + 10], 
              cex = .3, pch = 1, col = mcmc_nngp_list$Vecchia_approx_DAG$field_position$var_idx
        )
        
        #####
        # in order to check density : 
        # 1) set innovation of covparms to be 0
        # 2) make number of observations so small (~50) that there is only one cluster
        # 3) log diff should be 0 because it corresponds to 
        #    analytical sampling of the latent field
        print(
          ( 
            end_density-
              start_density
          )          +
            (
              reverse_field_proposal_density$transition_logdens
              -
                field_proposal$transition_logdens
            )
        )
        
        chain$params$field = field_proposal$field
        
        
        if(
          (
            # should be 0 if covariance parameters not moved
           ( 
             end_density-
            start_density
            )+
            (
              reverse_field_proposal_density$transition_logdens-
            field_proposal$transition_logdens
            )
          )>log(runif(1))
        ){
          print("tatato!")
          chain$params$field = field_proposal$field
          }
        plot(rep(stuff_for_plots$locs_no_na[,1], stuff_for_plots$n_var), stuff_for_plots$y_true[,,20], 
             col = rep(c("lightgray", "lightpink", "lightgreen", "lightblue", "lightcyan", "lightpink"), each = stuff_for_plots$n_loc), 
             cex = .3, pch = 15, 
             xlab = "spatial site", 
             ylab = "latent field and true field"
        )
        points(mcmc_nngp_list$locs[mcmc_nngp_list$Vecchia_approx_DAG$field_position$location_idx,1], 
               chain$params$field[,mcmc_nngp_list$useful_stuff$buffer_depth + 20], 
               cex = .3, pch = 1, col = mcmc_nngp_list$Vecchia_approx_DAG$field_position$var_idx
        )
        
      }
      # updating kernel
      chain$kernels$var_wise_ancillary[var_idx] = homeostasy(
        log_var = chain$kernels$var_wise_ancillary[var_idx], 
        accepted = accepted, target = .5, 
        2/sqrt(iter)
      )
    }

    
    if(iter/10==iter%/%10){
      par(mfrow = c(3,1))
      par(mar = c(2,2,2,2))
      plot(records$log_range_vec[1,,seq(match(0, records$log_range_vec[1,,])-1)], ylab = "")
      plot(records$log_range_vec[2,,seq(match(0, records$log_range_vec[1,,])-1)], ylab = "")
      plot(records$log_range_vec[3,,seq(match(0, records$log_range_vec[1,,])-1)], ylab = "")
      plot(records$smoothness_vec[1,,seq(match(0, records$smoothness_vec[1,,])-1)], ylab = "")
      plot(records$smoothness_vec[2,,seq(match(0, records$smoothness_vec[1,,])-1)], ylab = "")
      plot(records$smoothness_vec[3,,seq(match(0, records$smoothness_vec[1,,])-1)], ylab = "")
      par(mfrow = c(1,1))
    }
    ####################
    # Updating records #
    ####################
    for(name in names(chain$params)) records[[name]][,,iter] = chain$params[[name]]
  }
}


plot(records$smoothness_vec[1,,seq(match(0, records$log_range_vec[1,,])-1)], ylab = "smoothness 1")
plot(records$smoothness_vec[2,,seq(match(0, records$log_range_vec[1,,])-1)], ylab = "smoothness 2")
plot(records$smoothness_vec[3,,seq(match(0, records$log_range_vec[1,,])-1)], ylab = "smoothness 3")

plot(records$log_range_vec[1,, ])
plot(records$log_range_vec[2,, ])
plot(records$log_range_vec[3,, ])
plot(records$smoothness_vec[1,,])
plot(records$smoothness_vec[2,,])
plot(records$smoothness_vec[3,,])

dev.off()
par(mar = rep(2, 4))

ACFs = apply(records$field[,30,seq(match(0, records$log_range_vec[1,,])-1), drop = F], 1, function(x)c(acf(c(x), lag.max = 30, plot = F)$acf))
ACF_q = apply(ACFs, 1, function(x)quantile(x, c(0, .1, .25, .5, .75, .9, 1)))
tatato = acf(c(records$field[1,30,seq(match(0, records$log_range_vec[1,,])-1), drop = F]), plot = F)
plot(tatato, type = "n", ylim = c(-1, 1))
#apply(ACFs, 2, function(x)lines(seq(0, 30), x, col = "lightgray"))
apply(ACF_q, 1, function(x)lines(seq(0, 30), x))
