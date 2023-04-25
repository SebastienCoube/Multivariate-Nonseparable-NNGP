

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

get_noise_precisions = function(
    noise_info, 
    useful_stuff, 
    Vecchia_approx_DAG, 
    time_begin, 
    time_end
)
{
  # time periods
  k = seq(time_begin, time_end)
  noise_precision_matrices = 
    parallel::mclapply(
      mc.cores = parallel::detectCores()-1,
      X = k, 
      function(i)
      {   
        noise_precision_x = 
          unlist(lapply(seq(useful_stuff$n_loc), function(j)
          {
             noise_info[[j]][[i]]$tau_precision
          }))
        if(is.null(noise_precision_x))return(
          Matrix::sparseMatrix(
            i = 1,
            j = 1, 
            x = 0, 
            dims = c(useful_stuff$n_field, useful_stuff$n_field)
          )
        )
        Matrix::sparseMatrix(
          i = c(useful_stuff$noise_precision_i[[i]]),
          j = c(useful_stuff$noise_precision_j[[i]]), 
          x = c(noise_precision_x), 
          dims = c(useful_stuff$n_field, useful_stuff$n_field)
        )
      })
  
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
  Sys.time()-t1
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

get_B = function(locs, impulsion_centers, basis_function_range, Vecchia_approx_DAG){
  B = 
    fields::rdist(
      impulsion_centers, 
      locs[Vecchia_approx_DAG$field_position$location_idx,]
    )
  threshold = sqrt(-basis_function_range^2*log(0.01)/2)
  B[B> threshold]=0
  B[B!=0] = exp(-2*(basis_function_range)^-2*(B[B!=0])^2)
  B[B<.01]=0
  B = t(B)
  B = Matrix::sparseMatrix(
    i = rep(seq(nrow(B)), apply(B, 1, function(x)sum(x!=0))),
    j = unlist(useful_stuff$n_var_y * (unlist(apply(B, 1, function(x)which(x!=0)))-1) + 
                 unlist(mapply(function(x, y)rep(x, y), Vecchia_approx_DAG$field_position$var_idx, apply(B, 1, function(x)sum(x!=0))))), 
    x = c(unlist(apply(B, 1, function(x)x [x!=0]))), 
    dims = c(useful_stuff$n_field, useful_stuff$n_var_y*ncol(B)) 
  )
  B
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





