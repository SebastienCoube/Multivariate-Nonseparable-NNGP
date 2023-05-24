homeostasy  = function(log_var, accepted, target  = .25, eps = .01)
{
  if(accepted) log_var = log_var + eps * (1-target)
  if(!accepted) log_var = log_var - eps * target
  return(log_var)
}


do_100_updates = function(chain, mcmc_nngp_list, kernel_learning_rate, thinning_rate = 1)
{
  iter = 1
  records = lapply(chain$params, function(x)array(0, c(dim(x), 200)))
  
  for(iter in seq(iter, 200)){
    print(iter)
    #################################################################
    # Sampling the latent field using the "blocks 'n' bases" scheme #
    #################################################################
    if(iter == 1 | iter/2 ==iter %/% 2)
      t1 = Sys.time()
    {
      time_begin = mcmc_nngp_list$useful_stuff$time_depth
      time_end = mcmc_nngp_list$useful_stuff$n_time_periods - mcmc_nngp_list$useful_stuff$time_depth +1
      # fixing maximum cluster size
      cluster_size_target = 100
      # partitioning the space into clusters whose size is lower than cluster_size_target
      recursive_k_means_clust = recursive_k_means(locs_ = mcmc_nngp_list$useful_stuff$locs_repeated, cluster_size_target)
      basis_functions = c(
        # basis functions corresponding to usual 
        get_indicator_basis_functions(recursive_k_means_clust = recursive_k_means_clust, useful_stuff = mcmc_nngp_list$useful_stuff), 
        # groups of spatial basis functions
        unlist(lapply(
          get_grids(
            points = mcmc_nngp_list$useful_stuff$locs_repeated, 
            cluster_size_target = cluster_size_target), 
          function(x)get_basis_functions(
            points = mcmc_nngp_list$useful_stuff$locs_repeated, 
            tile_info = x, 
            cluster_size_target = cluster_size_target, 
            Vecchia_approx_DAG = mcmc_nngp_list$Vecchia_approx_DAG, 
            useful_stuff = mcmc_nngp_list$useful_stuff)),
          recursive = F)
      )
      # coloring to find parallelization groups
      coloring=  color_basis_functions(basis_functions, precision_blocks = chain$stuff$precision_blocks)
      colo = 1
      # looping on colors
      for(colo in seq(max(coloring)))
      { 
        #print(colo)
        
        # precision times latent field
        Q_field = vecchia_blocks_mult(
          x = chain$params$field,
          vecchia_blocks = chain$stuff$vecchia
        )
        Q_field[,-seq(time_begin, time_end)] = 0
        Q_field = vecchia_blocks_t_mult(
          x = Q_field,
          transposed_vecchia_blocks = chain$stuff$transposed_vecchia_blocks
        )
        Q_field = matrix(Q_field, nrow = mcmc_nngp_list$useful_stuff$n_field)
        # tau y - field
        y_minus_field = chain$params$field[,seq(time_begin, time_end)] - mcmc_nngp_list$useful_stuff$y_loc_var_format_no_NA[,seq(time_begin, time_end)]
        tau_y_minus_field = do.call(cbind, mapply(function(x,y) as.vector(Matrix::crossprod(x, y)), chain$stuff$noise_precisions[time_begin:time_end], split(y_minus_field, col(y_minus_field)), SIMPLIFY = F))
        i_basis_function=1
        eps = parallel::mclapply(
          mc.cores = 5, which(coloring == colo)
          ,
          function(i_basis_function)
          {
            # part with observations ####
            B_tau_y_minus_field = Matrix::crossprod(basis_functions[[i_basis_function]], tau_y_minus_field)
            B_tau_B_x = list()
            B_tau_B_i = list()
            B_tau_B_j = list()
            for(i_time in seq(time_begin, time_end)){
              B_tau_B = Matrix::crossprod(basis_functions[[i_basis_function]], Matrix::crossprod(chain$stuff$noise_precisions[time_begin:time_end][[i_time-time_begin+1]], basis_functions[[i_basis_function]]))
              B_tau_B = Matrix::forceSymmetric(B_tau_B)
              if(length(B_tau_B@x)>0)
              {
                B_tau_B_x = c(B_tau_B_x, list(B_tau_B@x))
                B_tau_B_i = c(B_tau_B_i, list((i_time-time_begin)*ncol(basis_functions[[i_basis_function]]) + B_tau_B@i+1))
                B_tau_B_j = c(B_tau_B_j, list((i_time-time_begin)*ncol(basis_functions[[i_basis_function]]) + findInterval(seq(length(B_tau_B@x))-1,B_tau_B@p[-1])+1))
              }
            }
            B_tau_B_x = unlist(B_tau_B_x)
            B_tau_B_i = unlist(B_tau_B_i)
            B_tau_B_j = unlist(B_tau_B_j)
            
            # part with prior ####
            B_Q_B = lapply(chain$stuff$precision_blocks, function(x)Matrix::crossprod(basis_functions[[i_basis_function]], x) %*% basis_functions[[i_basis_function]])
            B_Q_B[[1]]=Matrix::forceSymmetric(B_Q_B[[1]], uplo = "L")
            prior_precision = get_block_toeplitz(B_Q_B, time_end- time_begin+1)
            B_Q_field = (Matrix::crossprod(basis_functions[[i_basis_function]], Q_field[,seq(time_begin, time_end)]))
            
            # posterior ####
            posterior_chol = Matrix::chol(Matrix::sparseMatrix(
              i = c(prior_precision$i, B_tau_B_j), 
              j = c(prior_precision$j, B_tau_B_i), 
              x = c(prior_precision$x, B_tau_B_x), 
              symmetric = T, 
              dims = rep((time_end-time_begin+1)*ncol(basis_functions[[i_basis_function]]), 2)
            ))
            
            # sampling coefficients ####
            eps = 
              matrix(
                Matrix::solve(
                  posterior_chol, 
                  rnorm(length(B_tau_y_minus_field))
                  + Matrix::solve(
                    Matrix::t(posterior_chol), 
                    as.vector(
                      - B_Q_field 
                      - B_tau_y_minus_field
                    )
                  )
                ), 
                nrow = ncol(basis_functions[[i_basis_function]])
              )
            return(eps)
          }, mc.preschedule = F
        )
        
          chain$params$field[,time_begin:time_end] = 
            chain$params$field[,time_begin:time_end] +
            as.matrix(Reduce("+", mapply(function(e, b) b%*%e, eps, basis_functions[which(coloring==colo)])))
          
        #if(iter/10==iter%/%10)
          {
          plot(rep(stuff_for_plots$locs_no_na[,1], stuff_for_plots$n_var), stuff_for_plots$y_true[,,10], 
               col = rep(c("lightgray", "lightpink", "lightblue"), each = stuff_for_plots$n_loc), 
               cex = .3, pch = 15, 
               xlab = "spatial site", 
               ylab = "latent field and true field", 
               main = paste("pass", i, "color", colo)
          )
          points(mcmc_nngp_list$locs[mcmc_nngp_list$Vecchia_approx_DAG$field_position$location_idx,1], 
                 chain$params$field[,mcmc_nngp_list$useful_stuff$buffer_depth + 10], 
                 cex = .3, pch = 1, col = mcmc_nngp_list$Vecchia_approx_DAG$field_position$var_idx + (mcmc_nngp_list$Vecchia_approx_DAG$field_position$var_idx==3)
          )
        }
      }
    }
      Sys.time()-t1
    #########################
    # Covariance parameters REMEMBER TO ADD A #
    #########################
    var_idx = 1
    for(var_idx in seq(mcmc_nngp_list$useful_stuff$n_var_y))
    {
      print(paste("var_idx=", var_idx))
      #print(var_idx)
      #############
      # ancillary #
      #############
      proposal = #t(chol(matrix(c(1, -.8, -.8, 1), 2))) %*% 
        rnorm(2, sd = exp(.5*chain$kernels$var_wise_ancillary[var_idx]))
      #print(proposal)
      
      # proposing range
      proposed_log_range_vec = chain$params$log_range_vec
      proposed_log_range_vec[var_idx] = 
        proposed_log_range_vec[var_idx] +
        (
        mcmc_nngp_list$hierarchical_model$log_range_prior[var_idx, 2] - 
        mcmc_nngp_list$hierarchical_model$log_range_prior[var_idx, 1]
        )* proposal[1]
        #pseudo_sigmoid(
        #  pseudo_logit(
        #    proposed_log_range_vec[var_idx], 
        #    v_min = mcmc_nngp_list$hierarchical_model$log_range_prior[var_idx, 1],
        #    v_max = mcmc_nngp_list$hierarchical_model$log_range_prior[var_idx, 2]
        #  ) + proposal[1],
        #  v_min = mcmc_nngp_list$hierarchical_model$log_range_prior[var_idx, 1],
        #  v_max = mcmc_nngp_list$hierarchical_model$log_range_prior[var_idx, 2]
        #)
      
      # proposing smoothness
      proposed_smoothness_vec = chain$params$smoothness_vec
      proposed_smoothness_vec[var_idx] = 
        proposed_smoothness_vec[var_idx] +
        (
          mcmc_nngp_list$hierarchical_model$smoothness_prior[var_idx, 2] - 
            mcmc_nngp_list$hierarchical_model$smoothness_prior[var_idx, 1]
        )* proposal[2]
      ###   pseudo_sigmoid(
      ###     pseudo_logit(
      ###       proposed_smoothness_vec[var_idx],
      ###       v_min = mcmc_nngp_list$hierarchical_model$smoothness_prior[var_idx, 1],
      ###       v_max = mcmc_nngp_list$hierarchical_model$smoothness_prior[var_idx, 2]
      ###     ) + proposal[2] ,
      ###     v_min = mcmc_nngp_list$hierarchical_model$smoothness_prior[var_idx, 1],
      ###     v_max = mcmc_nngp_list$hierarchical_model$smoothness_prior[var_idx, 2]
      ###   )
      
      # checking prior
      accepted = F
      if(
        (proposed_log_range_vec[var_idx]>mcmc_nngp_list$hierarchical_model$log_range_prior[var_idx, 1])&
        (proposed_log_range_vec[var_idx]<mcmc_nngp_list$hierarchical_model$log_range_prior[var_idx, 2])&
        (proposed_smoothness_vec[var_idx]>mcmc_nngp_list$hierarchical_model$smoothness_prior[var_idx, 1])&
        (proposed_smoothness_vec[var_idx]<mcmc_nngp_list$hierarchical_model$smoothness_prior[var_idx, 2])
      ){
        #print("prior accepted")
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
        proposed_field = vecchia_blocks_solve(
          x = vecchia_blocks_mult(x = chain$params$field, vecchia_blocks = chain$stuff$vecchia_blocks), 
          vecchia_blocks = proposed_vecchia_blocks
        )
        
        
        
        
        if(iter/10==iter%/%10){
          plot(rep(stuff_for_plots$locs_no_na[,1], stuff_for_plots$n_var), stuff_for_plots$y_true[,,10], 
               col = rep(c("lightgray", "lightpink", "lightblue"), each = stuff_for_plots$n_loc), 
               cex = .3, pch = 15, 
               xlab = "spatial site", 
               ylab = "latent field and true field"
          )
          points(mcmc_nngp_list$locs[mcmc_nngp_list$Vecchia_approx_DAG$field_position$location_idx,1], 
                 chain$params$field[,mcmc_nngp_list$useful_stuff$buffer_depth + 10], 
                 cex = .3, pch = 1, col = mcmc_nngp_list$Vecchia_approx_DAG$field_position$var_idx + (mcmc_nngp_list$Vecchia_approx_DAG$field_position$var_idx==3)
          )
          
          plot(rep(stuff_for_plots$locs_no_na[,1], stuff_for_plots$n_var), stuff_for_plots$y_true[,,10], 
               col = rep(c("lightgray", "lightpink", "lightblue"), each = stuff_for_plots$n_loc), 
               cex = .3, pch = 15, 
               xlab = "spatial site", 
               ylab = "latent field and true field"
          )
          points(mcmc_nngp_list$locs[mcmc_nngp_list$Vecchia_approx_DAG$field_position$location_idx,1], 
                 proposed_field[,mcmc_nngp_list$useful_stuff$buffer_depth + 10], 
                 cex = .3, pch = 1, col = mcmc_nngp_list$Vecchia_approx_DAG$field_position$var_idx + (mcmc_nngp_list$Vecchia_approx_DAG$field_position$var_idx==3)
          )
        }
        
        if(
          -.5*(
            Reduce("+", 
                   mapply(
                     function(v, M)as.vector(t(v) %*% M %*% v),
                     v = split(mcmc_nngp_list$useful_stuff$y_loc_var_format_no_NA - proposed_field, col(proposed_field)),
                     M = chain$stuff$noise_precisions
                   )
            ) - 
              Reduce("+", 
                     mapply(
                       function(v, M)as.vector(t(v) %*% M %*% v),
                       v = split(mcmc_nngp_list$useful_stuff$y_loc_var_format_no_NA - chain$params$field, col(proposed_field)),
                       M = chain$stuff$noise_precisions
                     )
              ) 
          ) > log(runif(1))
        ){
          print("tatatooooo!")
          chain$params$field = proposed_field
          chain$params$log_range_vec = proposed_log_range_vec
          chain$params$smoothness_vec = proposed_smoothness_vec
          chain$stuff$vecchia_blocks = proposed_vecchia_blocks
          
          accepted = T
        }
      }
      # updating kernel
      chain$kernels$var_wise_ancillary[var_idx] = homeostasy(
        log_var = chain$kernels$var_wise_ancillary[var_idx], 
        accepted = accepted, target = .5, 
        2/sqrt(iter)
      )
      

      
      
      ##############
      # sufficient #
      ##############
      proposal =  
        rnorm(2, sd = exp(.5*chain$kernels$var_wise_sufficient[var_idx]))
      #print(proposal)
      
      # proposing range
      proposed_log_range_vec = chain$params$log_range_vec
      proposed_log_range_vec[var_idx] = 
        proposed_log_range_vec[var_idx] + 
        (
          mcmc_nngp_list$hierarchical_model$log_range_prior[var_idx, 2] - 
            mcmc_nngp_list$hierarchical_model$log_range_prior[var_idx, 1]
        )* proposal[1]
      #pseudo_sigmoid(
      #  pseudo_logit(
      #    proposed_log_range_vec[var_idx], 
      #    v_min = mcmc_nngp_list$hierarchical_model$log_range_prior[var_idx, 1],
      #    v_max = mcmc_nngp_list$hierarchical_model$log_range_prior[var_idx, 2]
      #  ) + proposal[1],
      #  v_min = mcmc_nngp_list$hierarchical_model$log_range_prior[var_idx, 1],
      #  v_max = mcmc_nngp_list$hierarchical_model$log_range_prior[var_idx, 2]
      #)
      
      # proposing smoothness
      proposed_smoothness_vec = chain$params$smoothness_vec
      proposed_smoothness_vec[var_idx] = 
        proposed_smoothness_vec[var_idx] + 
        (
          mcmc_nngp_list$hierarchical_model$smoothness_prior[var_idx, 2] - 
            mcmc_nngp_list$hierarchical_model$smoothness_prior[var_idx, 1]
        )* proposal[2]
      #####   pseudo_sigmoid(
      #####     pseudo_logit(
      #####       proposed_smoothness_vec[var_idx],
      #####       v_min = mcmc_nngp_list$hierarchical_model$smoothness_prior[var_idx, 1],
      #####       v_max = mcmc_nngp_list$hierarchical_model$smoothness_prior[var_idx, 2]
      #####     ) + proposal[2] ,
      #####     v_min = mcmc_nngp_list$hierarchical_model$smoothness_prior[var_idx, 1],
      #####     v_max = mcmc_nngp_list$hierarchical_model$smoothness_prior[var_idx, 2]
      #####   )
      
      
      # checking prior
      accepted = F
      if(
        (proposed_log_range_vec[var_idx]>mcmc_nngp_list$hierarchical_model$log_range_prior[var_idx, 1])&
        (proposed_log_range_vec[var_idx]<mcmc_nngp_list$hierarchical_model$log_range_prior[var_idx, 2])&
        (proposed_smoothness_vec[var_idx]>mcmc_nngp_list$hierarchical_model$smoothness_prior[var_idx, 1])&
        (proposed_smoothness_vec[var_idx]<mcmc_nngp_list$hierarchical_model$smoothness_prior[var_idx, 2])
      ){
        #print("prior accepted")
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
        
        if(
          (
            (mcmc_nngp_list$useful_stuff$n_time_periods - mcmc_nngp_list$useful_stuff$time_depth+1) * 
            (
              +sum(log(Matrix::diag(proposed_vecchia_blocks$triangular_on_diag)))
              -sum(log(Matrix::diag(chain$stuff$vecchia_blocks$triangular_on_diag)))
            )
            -.5*(
              sum(vecchia_blocks_mult(chain$params$field, proposed_vecchia_blocks)^2) - 
              sum(vecchia_blocks_mult(chain$params$field, chain$stuff$vecchia_blocks)^2)
            )
          )
          > log(runif(1))
          ){
          print("tatataaaaaaa!")
          chain$params$log_range_vec = proposed_log_range_vec
          chain$params$smoothness_vec = proposed_smoothness_vec
          chain$stuff$vecchia_blocks = proposed_vecchia_blocks
          
          accepted = T
        }
      }
      chain$kernels$var_wise_sufficient[var_idx] = homeostasy(
        log_var = chain$kernels$var_wise_sufficient[var_idx], 
        accepted = accepted, target = .5, 
        eps = 2/sqrt(iter)
      )
        ##################
        # Updating stuff #
        ##################
        chain$stuff$transposed_vecchia_blocks = vecchia_blocks_t(chain$stuff$vecchia_blocks)
        chain$stuff$precision_blocks = get_precision_blocks(vecchia_blocks = chain$stuff$vecchia_blocks)
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
