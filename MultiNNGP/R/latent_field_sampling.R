
B Q FIELD WRONG, CHECK CONDITIONAL FORMULA


# get precision blocks. Each block corresponds to a time lag, the first being lag = 0
# precision blocks starting from diag and going towards the bottom or the left

get_precision_blocks = function(vecchia_blocks)
{
  if(length(vecchia_blocks$rectangular_below_diag)==0) return(list(Matrix::crossprod(vecchia_blocks$triangular_on_diag)))
  res = list()
  res[[1]] = Matrix::crossprod(vecchia_blocks$triangular_on_diag) + Reduce("+", lapply(vecchia_blocks$rectangular_below_diag, Matrix::crossprod))
  max_lag = length(vecchia_blocks$rectangular_below_diag)
  for(i_lag in seq(max_lag))
  {
    res[[i_lag+1]] = Matrix::crossprod(vecchia_blocks$triangular_on_diag, vecchia_blocks$rectangular_below_diag[[max_lag-i_lag+1]])
    if(max_lag!=i_lag) 
      res[[i_lag+1]] = res[[i_lag+1]] +  
        Reduce("+", 
               mapply(
                 Matrix::crossprod, 
                 vecchia_blocks$rectangular_below_diag[seq(i_lag+1, max_lag)],
                 vecchia_blocks$rectangular_below_diag[seq(max_lag-i_lag) ]
                 )
               )
  }
  res
}

get_block_toeplitz = function(blocks, block_dim)
{
  blocks= lapply(blocks, function(x)as(x, "TsparseMatrix"))
  n = ncol(blocks[[1]])
  i_idx = 
    unlist(
    lapply(
      seq(length(blocks)),
      function(i)
      {
        unlist(lapply(seq(i, block_dim) - 1, function(j)blocks[[i]]@i+1 + j*n))
      }
    )
  )  
  j_idx = 
    unlist(
    lapply(
      seq(length(blocks)),
      function(i){
        unlist(lapply(seq(block_dim - i+1) - 1, function(j)blocks[[i]]@j+1 + j*n ))
      }
    )
  )  
  x = unlist(
    mapply(
      function(i, x)
      {
        rep(x, block_dim-i+1)
      }, 
      i = seq(length(blocks)), 
      x = lapply(blocks, function(M)(M@x))
    )
  ) 
  return(list("i" = i_idx, "j" = j_idx, "x" = x))
}

get_precision_block_chol = function(B_tau_B, B_Q_B)
{
  # Recursive use of :
  #            |                      |          
  #            |                      |          
  #     A      | B            chol(A) | chol(A)-t B          
  #            |     ---->            |             
  #   ---------|---          ---------|---          
  #     Bt     | C               0    | chol(C-Bt chol(A)-1 B) 
  #   
  # For chol(A)-t B :  Here we have 
  #                         |               | 0 |    
  #                we don't | 0             | 0 |    
  #  chol(A)-t B =    care  |               | 0 |    
  #                ---------|------         |---|    
  #                we don't | chol(A22)-t   | B2|       
  #                   care  |               | B2|  
  #                                              
  time_depth = length(B_Q_B)
  res = list()
  i_time = 1
  transposed_precisions = lapply(B_Q_B, Matrix::t)
  for(i_time in seq(length(B_tau_B)))
  {
    print(i_time)
    res[[i_time]] = list()
    # chol(A)-t B
    diag_block = B_Q_B[[1]] + B_tau_B[[i_time]]
    
    if(i_time>1)
    {
      condition_depth = min(time_depth-1, i_time-1)
      # starting by last block
      res[[i_time]]$blocks_above_diag = transposed_precisions[seq(2, 1+condition_depth)]
      for(i_lag in seq(condition_depth, 1))
      {
        #i_lag_reverse_order = condition_depth - i_lag + 1
        # solving by triangular on diag
        res[[i_time]]$blocks_above_diag[[i_lag]] = backsolve(
          r = res[[i_time-i_lag]]$triangular_on_diag, 
          x = res[[i_time]]$blocks_above_diag[[i_lag]], transpose = T)
        # updating following blocks
        if(i_lag>1)
        {
          for(i_lag_ in seq(1, i_lag-1))
          {
            res[[i_time]]$blocks_above_diag[[i_lag_]] = 
              res[[i_time]]$blocks_above_diag[[i_lag_]] - 
              t(res[[i_time-i_lag_]]$blocks_above_diag[[i_lag-i_lag_]]) %*%
              res[[i_time]]$blocks_above_diag[[i_lag]]
          }
        }
      }
      names(res[[i_time]]$blocks_above_diag) = paste("lag=", seq(condition_depth), sep = "")
      # removing Bt chol(A)-1 B
      diag_block = diag_block - Reduce("+", lapply(res[[i_time]]$blocks_above_diag, crossprod))
    }
    res[[i_time]]$triangular_on_diag =chol(diag_block)
  }
}


DSATUR = function(M)
{
  n = ncol(M)
  #deducting degrees
  degrees = as.vector(rep(1, n)%*%M)
  #getting adjacent nodes of a given node
  neighbors = split(M@i+1, rep(seq_along(diff(M@p)),diff(M@p)))
  #creating a color * node matrix of incompatibilities
  incompatibilities = matrix(0, nrow(M)+1, max(degrees))
  # creating color vector
  cols = rep(0, nrow(M))
  for(i in seq(nrow(M)))
  {
    idx = which(cols == 0)
    # computing saturation degrees
    dsatur = sapply(idx, function(idx)length(setdiff(unique(cols[neighbors[[idx]]]), 0)))
    # finding maximum saturation
    max_satur = which(dsatur == max(dsatur))
    # finding maximum degree among max satured
    max_degree = which.max(degrees[idx][max_satur])
    # find new colored location -1 max satur -2 max degree
    colored_location =  idx[max_satur][max_degree][1]
    cols[colored_location] = match(0, incompatibilities[colored_location,])
    incompatibilities[neighbors[[colored_location]],cols[colored_location]] = 1
  }
  return(cols)  
}


naive_greedy_coloring = function(M)
{
  M[2,1] = M[2,1]+1
  M[2,1] = M[2,1]+-1
  #number of nodes
  n_obs = nrow(M)
  #deducting degrees
  degrees = as.vector(rep(1, n_obs)%*%M)
  #getting adjacent nodes of a given node
  idx = split(M@i+1, rep(seq_along(diff(M@p)),diff(M@p)))
  #creating a color * node matrix of incompatibilities
  incompatibilities = matrix(0, n_obs+1, max(degrees))
  
  cols = rep(0, n_obs)
  
  for(i in seq(n_obs))
  {
    cols[i] = match(0, incompatibilities[i,])
    incompatibilities[idx[[i]],cols[i]] = 1
  }
  return(cols)
}

## field = chain$params$field
## epsilon = NULL
## precision_blocks = chain$stuff$precision_blocks
## transposed_vecchia_blocks = chain$stuff$transposed_vecchia_blocks
## vecchia_blocks = chain$stuff$vecchia_blocks
## noise_precisions = chain$stuff$noise_precisions
## 
## 
## list2env(
##   get_blocks_n_bases_with_witnesses(
##     mcmc_nngp_list = mcmc_nngp_list, 
##     space_cluster_size_target = 120, 
##     chain = chain),
##   envir =  environment()
## )
## list2env(
##   get_time_split(
##     time_depth = mcmc_nngp_list$useful_stuff$time_depth, 
##     n_time_periods = mcmc_nngp_list$useful_stuff$n_time_periods, 
##     time_target_size = 100),
##   envir =  environment()
## )
## 
## remove(chain)

# either proposes a new latent field or
# evaluates transition density between two latent fields
latent_field_blocked_sampling = function(
    field, 
    epsilon = NULL, # if NULL then new latent field is sampled, 
    # else reverse transition ratio is computed
    mcmc_nngp_list, 
    time_coloring, spatial_coloring, 
    time_split, basis_functions,
    precision_blocks, 
    vecchia_blocks, 
    transposed_vecchia_blocks, 
    noise_precisions
)
{
  # reversing epsilon 
  if(!is.null(epsilon))epsilon = lapply(epsilon, function(x)lapply(x, function(x)-x))
  # creating epsilon
  if(is.null(epsilon)){
    epsilon = 
      lapply(seq(nrow(time_split)), 
             function(x)lapply(basis_functions, function(x)NULL))
  }
  # initializing transition log density
  transition_logdens = 0
  # for testing 
  time_color = 1
  spatial_color = 1
  # looping over time and space coloring
  for(time_color in seq(max(time_coloring)))
  {
    for(spatial_color in seq(max(spatial_coloring)))
    { 
      # contribution of observations to field mean
      y_minus_field = field - mcmc_nngp_list$useful_stuff$y_loc_var_format_no_NA
      tau_y_minus_field = do.call(cbind, mapply(function(x,y) as.vector(Matrix::crossprod(x, y)), noise_precisions, split(y_minus_field, col(y_minus_field)), SIMPLIFY = F))
      # index of selected spatial basis functions and time periods
      space_time_idx = unlist(
        recursive = F,
        lapply(which(spatial_coloring==spatial_color), 
               function(x)lapply(which(time_coloring == time_color), function(y)c(x,y))
        )
      )
      space_time_couple = space_time_idx[[1]]
      transitions = parallel::mclapply(
        mc.cores = 8, space_time_idx
        ,
        function(space_time_couple)
        {
          i_basis_function = space_time_couple[[1]]
          i_time_split = space_time_couple[[2]]
          time_begin = time_split[i_time_split,1]
          time_end   = time_split[i_time_split,2]
          
          # part with observations ####
          # contribution to the mean 
          B_tau_y_minus_field = Matrix::crossprod(basis_functions[[i_basis_function]], tau_y_minus_field[,seq(time_begin, time_end)])
          # contribution to precision 
          B_tau_B_x = list()
          B_tau_B_i = list()
          B_tau_B_j = list()
          for(i_time in seq(time_begin, time_end)){
            B_tau_B = Matrix::crossprod(basis_functions[[i_basis_function]], Matrix::crossprod(noise_precisions[[i_time]], basis_functions[[i_basis_function]]))
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
          # contribution to precision
          B_Q_B = lapply(precision_blocks, function(x)Matrix::crossprod(basis_functions[[i_basis_function]], x) %*% basis_functions[[i_basis_function]])
          B_Q_B[[1]]=Matrix::forceSymmetric(B_Q_B[[1]], uplo = "L")
          prior_precision = get_block_toeplitz(B_Q_B, time_end- time_begin+1)
          
          # contribution to the mean
          # precision times latent field
          Q_field = vecchia_blocks_mult(
            x = field[,
                      seq(
                        time_begin-mcmc_nngp_list$useful_stuff$time_depth+1, 
                        time_end+mcmc_nngp_list$useful_stuff$time_depth-1)], # time markov blanket
            vecchia_blocks = vecchia_blocks
          )
          Q_field = Q_field[,
                            -c(
                              seq(mcmc_nngp_list$useful_stuff$time_depth-1), 
                              seq((ncol(Q_field)- mcmc_nngp_list$useful_stuff$time_depth+2), ncol(Q_field)))]
          # subsetting in order to keep only places where basis functions are non null
          row_idx = which(as.vector(basis_functions[[i_basis_function]]%*% rep(1, ncol(basis_functions[[i_basis_function]]))!=0))
          B_transposed_vecchia_blocks = transposed_vecchia_blocks
          B_transposed_vecchia_blocks$t_triangular_on_diag = 
            Matrix::crossprod(basis_functions[[i_basis_function]],
                              B_transposed_vecchia_blocks$t_triangular_on_diag)
          B_transposed_vecchia_blocks$t_rectangular_below_diag = 
            lapply(B_transposed_vecchia_blocks$t_rectangular_below_diag, 
                   function(x)Matrix::crossprod(basis_functions[[i_basis_function]],x))
          B_Q_field = vecchia_blocks_t_mult(
            x = Q_field,
            transposed_vecchia_blocks = B_transposed_vecchia_blocks
          )
          
          # posterior #####
          
          t1 = Sys.time()
          posterior_chol = Matrix::chol(
            Matrix::sparseMatrix(
              i = c(prior_precision$i, B_tau_B_j), 
              j = c(prior_precision$j, B_tau_B_i), 
              x = c(prior_precision$x, B_tau_B_x), 
              symmetric = T, 
              dims = rep((time_end-time_begin+1)*ncol(basis_functions[[i_basis_function]]), 2)
            )
          )
          Sys.time()-t1
          # needs to be solved by cond_mean_precursor to give conditional mean
          cond_mean_precursor = Matrix::solve(
            Matrix::t(posterior_chol), 
            as.vector(
              - B_Q_field 
              - B_tau_y_minus_field
            )
          )
          # sampling coefficients if field is sampled ####
          eps = NULL
          if(is.null(epsilon[[i_time_split]][[i_basis_function]])){
            eta = 
0* 
              rnorm(length(cond_mean_precursor))
            
            eps =
              matrix(
                Matrix::solve(
                  posterior_chol, 
                  eta
                  + cond_mean_precursor
                ), 
                nrow = ncol(basis_functions[[i_basis_function]])
              )
          }
          # evaluating backward density if field is not sampled ####
          if(!is.null(epsilon[[i_time_split]][[i_basis_function]])){
            # solving equation : eps = solve(posterior_chol, eta + cond_mean_precursor)
            eta = posterior_chol %*% c(epsilon[[i_time_split]][[i_basis_function]]) - cond_mean_precursor
          }
          # transition density #### 
          transition_density = sum(log(Matrix::diag(posterior_chol)))  - # determinant of posterior precision
            .5 * sum(eta^2) 
          # returning
          return(
            list("eps" = eps, "transition_density" = transition_density)
          )
        }, mc.preschedule = F
      )
      
      #t1 = Sys.time()
      eps_idx = 1
      for(eps_idx in seq(length(transitions)))
      {
        transition_logdens = transition_logdens + transitions[[eps_idx]]$transition_density
        i_time_split = space_time_idx[[eps_idx]][2]
        i_basis_function = space_time_idx[[eps_idx]][1]
        if(is.null(epsilon[[i_time_split]][[i_basis_function]]))
        {
          epsilon[[i_time_split]][[i_basis_function]] = 
            transitions[[eps_idx]]$eps
        }
        time_begin = time_split[i_time_split,1]
        time_end   = time_split[i_time_split,2]
        field[,time_begin:time_end] = 
          as.matrix(
            field[,time_begin:time_end] +
              basis_functions[[i_basis_function]] %*% epsilon[[i_time_split]][[i_basis_function]]
          )
      }
      ##print(Sys.time()-t1)
      ## plot(rep(stuff_for_plots$locs_no_na[,1], stuff_for_plots$n_var), stuff_for_plots$y_true[,,80], 
      ##      col = rep(c("lightgray", "lightpink", "lightgreen", "lightblue", "lightcyan", "lightpink"), each = stuff_for_plots$n_loc), 
      ##      cex = .3, pch = 15, 
      ##      xlab = "spatial site", 
      ##      ylab = "latent field and true field"
      ## )
      ## points(mcmc_nngp_list$locs[mcmc_nngp_list$Vecchia_approx_DAG$field_position$location_idx,1], 
      ##        field[,mcmc_nngp_list$useful_stuff$buffer_depth + 80], 
      ##        cex = .3, pch = 1, col = mcmc_nngp_list$Vecchia_approx_DAG$field_position$var_idx
      ## )
    }
  }
  return(list("field" = field, "epsilon" = epsilon, "transition_logdens" = transition_logdens))
}


covparms_and_latent_field_density = 
  function(field, 
           vecchia_blocks, 
           noise_precisions, 
           mcmc_nngp_list
           ){
    y_minus_field = field - mcmc_nngp_list$useful_stuff$y_loc_var_format_no_NA
    return(
      (mcmc_nngp_list$useful_stuff$n_time_periods - mcmc_nngp_list$useful_stuff$time_depth+1) * sum(log(Matrix::diag(vecchia_blocks$triangular_on_diag)))
      -.5 * sum(vecchia_blocks_mult(x = field, vecchia_blocks = vecchia_blocks)^2) 
      -.5 * Reduce("+", mapply(function(x,y) as.vector(y %*% Matrix::crossprod(x, y)), noise_precisions, split(y_minus_field, col(y_minus_field)), SIMPLIFY = F))
    )
  }

