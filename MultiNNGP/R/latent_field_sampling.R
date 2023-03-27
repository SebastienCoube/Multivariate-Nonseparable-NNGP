

# get precision blocks. Each block corresponds to a time lag, the first being lag = 0
# precision blocks starting from diag and going towards the bottom or the left

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
get_block_toeplitz = function(blocks, block_dim)
{
  n_blocks = length(blocks)
  k = ncol(blocks[[1]])
  i_idx = list()
  j_idx = list()
  x = list()
  i_idx[[1]] = outer(row(blocks[[1]])[lower.tri(blocks[[1]], diag = T)], seq(0, block_dim-1)*k, "+")
  j_idx[[1]] = outer(col(blocks[[1]])[lower.tri(blocks[[1]], diag = T)], seq(0, block_dim-1)*k, "+")
  x[[1]]     = rep  (    blocks[[1]] [lower.tri(blocks[[1]], diag = T)], block_dim)
  
  if(block_dim==1)return(list("i" = i_idx, "j" = j_idx, "x" = x))
  for(i_block in seq(2, n_blocks))
  {
    i_idx[[i_block]] = (i_block-1)*k + outer(row(blocks[[i_block]]), seq(0, block_dim-i_block)*k, "+")
    j_idx[[i_block]] =                 outer(col(blocks[[i_block]]), seq(0, block_dim-i_block)*k, "+")
    x    [[i_block]] =                 rep  (    blocks[[i_block]],  block_dim-i_block+1)
  }
  i_idx = unlist(i_idx)
  j_idx = unlist(j_idx)
  x = unlist(x)
  return(list("i" = i_idx, "j" = j_idx, "x" = x))
}

##tatato = Matrix::sparseMatrix(
##  i = i_idx,
##  j = j_idx,
##  x = x, symmetric = T
##)


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


















useful_stuff$position_in_y_at_least_one_obs
get_BT_tau = function(noise_info, t_start, t_end, useful_stuff, Vecchia_approx_DAG, selected_locs)
{
  
  B_tau = lapply(
    seq(t_end-t_start), function(i)
    {
      matrix(0, nrow = useful_stuff$n_field, ncol = useful_stuff$n_var_y)
    }
  )
  for(i_time in seq(t_end-t_start+1)){
    for(i_site in selected_locs){
      i_field = Vecchia_approx_DAG$field_position$loc_match[[i_site]]
      i_var = Vecchia_approx_DAG$field_position$var_idx[i_field]
      if(length(useful_stuff$position_in_y_at_least_one_obs[[i_time]][[i_site]])>0)
      {
        B_tau[[i_time]][i_field, i_var[useful_stuff$position_in_y_at_least_one_obs[[i_time]][[i_site]]]] = 
          noise_info[[i_site]][[i_time]]$tau_precision
      }
    }
  }
  
  x = unname(unlist(lapply(seq(t_start, t_end), function(i_time)
  {
    lapply(noise_info[selected_locs], function(x)x[[i_time]]$tau_precision)
  })))
  j = unname(unlist(lapply(seq(t_start, t_end), function(i_time)
  {
    lapply(Vecchia_approx_DAG$field_position$loc_match[selected_locs], function(i_field)rep(Vecchia_approx_DAG$field_position$var_idx[i_field]))
  })))
    i_field = Vecchia_approx_DAG$field_position$loc_match[[i_site]]
    i_var = Vecchia_approx_DAG$field_position$var_idx[i_field]
  
  useful_stuff$position_in_y_at_least_one_obs[t_start:t_end]
  # precisions at selected locs and times
  lapply(noise_info[selected_locs], function(x)lapply(x[t_start:t_end], function(x)x$tau_precision))
}


get_BT_tau = function(noise_info, t_start, t_end, useful_stuff, Vecchia_approx_DAG, selected_locs)
{
  t1 = Sys.time()
  x = unname(unlist(
    lapply(
      noise_info[selected_locs], # selecting sites
      function(x)
        lapply(
          x[seq(t_start, t_end)], # selecting times
          function(x)x$tau_precision # selecting precision
        )
    )
  ))
  Sys.time()-t1
  t1 = Sys.time()
  unname(unlist(
    lapply(selected_locs, function(i_loc)
      lapply(seq(t_start, t_end), function(i_time)
      {
        i_idx = 
          seq(useful_stuff$n_var_y)[ # all possible variables of y
            useful_stuff$y_at_least_one_obs[[i_loc]][ # selecting variables of the field at site i_loc
              useful_stuff$position_in_y_at_least_one_obs[[i_time]][[i_loc]] # selecting variables that are observed at i_loc, i_time
            ]
          ]
        return(NULL)
        outer(i_idx, rep(1, length(i_idx)))
      }
      )
    )
  ))
  Sys.time()-t1
}

get_tau_y = function(y, useful_stuff, noise_info)
{
  tatato =lapply(noise_info[[1]], function(x)x$tau_precision)
  tatato = tatato[-which(sapply(tatato, is.null))]
  Matrix::bdiag(tatato)
}




vecchia_blocks = mcmc_nngp_list$chains[[1]]$stuff$vecchia
precision_blocks = get_precision_blocks(vecchia_blocks)
noise_info = mcmc_nngp_list$chains[[1]]$stuff$noise_info
X_noise_list = mcmc_nngp_list$covariates$X_noise
diag_precision_block = precision_blocks[[1]]
useful_stuff = mcmc_nngp_list$useful_stuff

loc_subset_idx = seq(useful_stuff$n_loc/2)


# evaluates 

get_field_obs_density = function(selected_sites, time_window, 
                                 field_subset, y_, noise_info, useful_stuff,
                                 Vecchia_approx_DAG
){
  -.5*sum(
    sapply(selected_sites, function(i_site)
      sum(
        sapply(seq(time_window[1], time_window[2]), function(i_time)
        {
          if(all(is.na(y_[i_site,,i_time])))return(0)
          idx =  which(!is.na(y_[i_site,,i_time]))
          return(
            sum(
              (
                y_[i_site,,i_time][idx] 
                - field_subset[Vecchia_approx_DAG$field_position$loc_match[[i_site]], i_time - time_window[1]+1][useful_stuff$position_in_y_at_least_one_obs[[i_time]][[i_site]]]
              ) *
                (
                  noise_info[[i_site]][[i_time]]$tau_precision %*%
                    (
                      y_[i_site,,i_time][idx] 
                      - field_subset[Vecchia_approx_DAG$field_position$loc_match[[i_site]], i_time - time_window[1]+1][useful_stuff$position_in_y_at_least_one_obs[[i_time]][[i_site]]]
                    ) 
                )
            )
          )
        }
        )
      )
    )
  )
}



## get_loc_var_pairs_partition = function(locs, Vecchia_approx_DAG, clust_size = 500)
## {
##   clust = rep(1, nrow(locs))
##   i=1
##   tablita = table(rep(clust, unlist(useful_stuff$n_field_per_site)))
##   while(any(tablita>clust_size))
##   {
##     for(cl in names(tablita)[tablita>clust_size])
##     {
##       n_clust =floor(5+runif(1, 0, 5))
##       if(n_in_clust<clust_size*5) n_clust = 2+rbinom(1,1,.5)
##       if(n_in_clust<clust_size*2) n_clust = 2
##       clust[clust == cl] = paste(i, kmeans(locs[clust==cl,], n_clust)$cl)
##       i=i+1
##     }
##     tablita = table(rep(clust, unlist(useful_stuff$n_field_per_site)))
##   }
##   clust = as.integer(as.factor(clust))
##   clust
## }
## 
## loc_var_pairs_partition = get_loc_var_pairs_partition(locs = mcmc_nngp_list$locs, Vecchia_approx_DAG = Vecchia_approx_DAG, clust_size = 400)
## plot(mcmc_nngp_list$locs[mcmc_nngp_list$Vecchia_approx_DAG$field_position$location_idx,], col = loc_var_pairs_partition)
## hist(table(loc_var_pairs_partition))
##  
##  
## ##get_clust_Markov_graph = function(Vecchia_approx_DAG, locs_partition)
## ##{
## ## field_dag_mat = Matrix::sparseMatrix(
## ##   mapply(c, 
## ##          Vecchia_approx_DAG$DAG$parents_same_time, 
## ##          Vecchia_approx_DAG$DAG$parents_previous_times 
## ##   )
## ##   )
## ##}
##  
##  
## selected_loc_var_pairs = which(loc_var_pairs_partition==3)
## sampling_time_depth = 8
## 
## get_prior_precision_cholesky = function(
##     precision_blocks, 
##     X_noise_list, 
##     useful_stuff, 
##     selected_loc_var_pairs, 
##     sampling_time_depth = 8, 
##     Vecchia_approx_DAG
## )
## {
##   
##   n_selected_field = length(selected_loc_var_pairs)
##   precision_blocks_ = lapply(precision_blocks, function(x)x[selected_loc_var_pairs,selected_loc_var_pairs])
##   below_diag_blocks = precision_blocks_; below_diag_blocks[[1]]=NULL
##   diag_precision_block = precision_blocks_[[1]]
##   
##   
##   block_diag_ij = sparseHessianFD::Matrix.to.Coord(diag_precision_block)
##   block_diag_i = block_diag_ij$rows[block_diag_ij$rows>=block_diag_ij$cols]
##   block_diag_j = block_diag_ij$cols[block_diag_ij$rows>=block_diag_ij$cols]
##   block_diag_x = diag_precision_block[cbind(block_diag_i, block_diag_j)]
##   block_diag_i = unlist(lapply(seq(0, sampling_time_depth-1), function(i)i*n_selected_field + block_diag_i))
##   block_diag_j = unlist(lapply(seq(0, sampling_time_depth-1), function(i)i*n_selected_field + block_diag_j))
##   block_diag_x = rep(block_diag_x, sampling_time_depth)
##   
##   # indices of below diag
##   below_diag_ij = lapply(below_diag_blocks, sparseHessianFD::Matrix.to.Coord)
##   below_diag_i = lapply(below_diag_ij, function(x)x$row)
##   below_diag_j = lapply(below_diag_ij, function(x)x$col)
##   below_diag_x = mapply(
##     function(mat, i, j)mat[cbind(i, j)],
##     below_diag_blocks,
##     below_diag_i, 
##     below_diag_j, 
##     SIMPLIFY = F
##   )
##   
##  M = matrix(0, sampling_time_depth, sampling_time_depth)
##  selected_blocks = (row(M)-col(M)>0)&(row(M)-col(M)<useful_stuff$time_depth)
##  selected_lags = (row(M)-col(M))[selected_blocks]
##  below_block_diag_x = unlist(below_diag_x[selected_lags])
##  below_block_diag_i = unlist(mapply(
##    function(lag, time_start)
##    {
##      below_diag_i[[lag]] + time_start * n_selected_field
##    },
##    selected_lags, 
##    (row(M)-1)[selected_blocks]
##  ))
##  below_block_diag_j = unlist(mapply(
##    function(lag, time_start)
##    {
##      below_diag_j[[lag]] + time_start * n_selected_field
##    },
##    selected_lags, 
##    (col(M)-1)[selected_blocks]
##  ))
##  precision_subset =  Matrix::sparseMatrix(
##    i = c(block_diag_i, below_block_diag_i),
##    j = c(block_diag_j, below_block_diag_j),
##    x = c(block_diag_x, below_block_diag_x), 
##    symmetric = T
##  )
##  Matrix::Cholesky(
##     precision_subset
##   )
## }
## 
## 
## get_prior_precision_choleskys = function(
##     precision_blocks, 
##     X_noise_list, 
##     useful_stuff, 
##     loc_var_pairs_partition, 
##     sampling_time_depth, 
##     Vecchia_approx_DAG
##   )
## {
##   parallel::mclapply(mc.cores = 10,
##     unique(loc_var_pairs_partition), function(i)
##     get_prior_precision_cholesky(
##       precision_blocks, 
##       X_noise_list, 
##       useful_stuff, 
##       selected_loc_var_pairs = which(loc_var_pairs_partition==i), 
##       sampling_time_depth, 
##       Vecchia_approx_DAG
##     )
##     )
## }
## 
## prior_precision_choleskys = 
## get_prior_precision_choleskys (
##     precision_blocks, 
##     X_noise_list, 
##     useful_stuff, 
##     loc_var_pairs_partition, 
##     sampling_time_depth, 
##     Vecchia_approx_DAG
##   )
## 
## useful_stuff$n_time_periods
## useful_stuff$time_depth
## sampling_time_depth
## time_begin = 30
## time_end = 38
## 
## field = mcmc_nngp_list$chains$chain_1$params$field
## 
## 
## precision_blocks_t = lapply(precision_blocks, Matrix::t)
## field = mcmc_nngp_list$chains$chain_1$params$field
## get_QAB_XB_field = function(
##   precision_blocks,   
##   precision_blocks_t,   
##   useful_stuff, 
##   field, 
##   selected_loc_var_pairs, 
##   time_begin, 
##   time_end
##   ){
##   interest_time_length = time_end-time_begin+1
##   field[time_begin:time_end, selected_loc_var_pairs]=0
##   field = field[seq(time_begin-useful_stuff$time_depth+1, time_end+useful_stuff$time_depth-1),]
##   res = matrix(0, interest_time_length, length(selected_loc_var_pairs))
##   for(i_lag in seq(-useful_stuff$time_depth+1, +useful_stuff$time_depth-1))
##   {
##     if(i_lag<=0){
##       res = res + as.matrix( field[seq(useful_stuff$time_depth + i_lag, interest_time_length + useful_stuff$time_depth - 1 + i_lag),]%*% precision_blocks[[abs(i_lag)+1]][,selected_loc_var_pairs])
##     }
##     if(i_lag>0){
##       res = res + as.matrix( field[seq(useful_stuff$time_depth + i_lag, interest_time_length + useful_stuff$time_depth - 1 + i_lag),]%*% precision_blocks_t[[abs(i_lag)+1]][,selected_loc_var_pairs])
##     }
##   }
##   res
## }
## 
## 
## # part of the precision matrix off the diagonal, between the selected field and the observations corresponding to them
## # In virtue of the data model, the rest of the observations is independent on the considered latent variables CONDITIONALLY ON the rest of the latent variables
## get_field_obs_precision = function(
##     noise_info, 
##     X_noise_list, 
##     useful_stuff, 
##     selected_loc_var_pairs, 
##     selected_times, 
##     Vecchia_approx_DAG,
##     y
## ){
##   
##   
##   selected_locs = unique(Vecchia_approx_DAG$field_position$location_idx[selected_loc_var_pairs])
##   selected_vars_per_loc = split(Vecchia_approx_DAG$field_position$var_idx[selected_loc_var_pairs], Vecchia_approx_DAG$field_position$location_idx[selected_loc_var_pairs])
##   useful_stuff$y_at_least_one_obs
##   noise_info
##   
##   #loc-var pairs of the DAG selected by loc only partition
##   n_field_per_site_subset = unlist(useful_stuff$n_field_per_site[loc_subset_idx])
##   cumsum_n_field_per_site_subset = cumsum(c(0, n_field_per_site_subset))
##   # indices of additional precision due to observations
##   additional_precision_indices = do.call(rbind, lapply(seq(length(loc_subset_idx)), function(j)
##   {
##     if(n_field_per_site_subset[j]==0)return(matrix(0, 0, 2))
##     
##     cbind(
##       cumsum_n_field_per_site_subset[j] + c(outer(seq(   n_field_per_site_subset[j]), rep(1, dim(y)[2]))), 
##       dim(y)[2]*(j-1) + c(outer(rep(1, n_field_per_site_subset[j]), seq(   dim(y)[2]))) 
##     )
##   }
##   )
##   )
##   additional_precision_i = additional_precision_indices[,1]
##   additional_precision_j = additional_precision_indices[,2]
##   # time periods
##   k = 1
##   if(useful_stuff$n_time_periods>1) k = seq(useful_stuff$buffer_depth+1, useful_stuff$n_time_periods)
##   field_obs_precision = 
##     parallel::mclapply(
##       mc.cores = parallel::detectCores()-1,
##       X = k, 
##       function(i)
##       {   
##         additional_precision_x = 
##           unlist(lapply(seq(length(loc_subset_idx)), function(j)
##           {
##             print(j)
##             if(n_field_per_site_subset[j]==0)return(NULL)
##             # creating matrix of 0 whose size corresponds to the number of latent field occurences at the selected location
##             # filling the matrix with noise precision, at the indices corresponding to observed variables
##             # eg: variables 1, 2, 4 are observed at site s. So tau_precision will have size 3*3. 
##             # But at time t, only variables 1 and 4 are observed. So tau_precision will look like: 
##             #    ?  0  ?
##             #    0  0  0
##             #    ?  0  ?
##             tau_precision = matrix(0, n_field_per_site_subset[j], dim(y)[2])
##             ##if(length(useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]])>1)
##             #{
##             #tau_precision[,useful_stuff$y_at_least_one_obs[[loc_subset_idx[j]]], drop = F][useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]],useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]]] = c(noise_info[[loc_subset_idx[[j]]]][[i]]$tau_precision)
##             tau_precision[useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]] , seq(ncol(tau_precision))[useful_stuff$y_at_least_one_obs[[loc_subset_idx[j]]]] [useful_stuff$position_in_y_at_least_one_obs[[i]][[loc_subset_idx[j]]]] ] = c(noise_info[[loc_subset_idx[[j]]]][[i]]$tau_precision)
##             #}
##             
##             c(tau_precision)
##           }))
##         Matrix::sparseMatrix(
##           i = c(additional_precision_i),
##           j = c(additional_precision_j), 
##           x = c(additional_precision_x)
##         ) 
##       })
## }
## 
## 
## 
## get_field_obs_precisions = function(
##     noise_info, 
##     X_noise_list, 
##     useful_stuff, 
##     loc_var_pairs_partition, 
##     Vecchia_approx_DAG,
##     y
## ){
##   #parallel::mc
##   lapply(#mc.cores = 10,
##                      unique(loc_var_pairs_partition), function(i)
##                        get_prior_precision_cholesky(
##                          noise_info, 
##                          X_noise_list, 
##                          useful_stuff, 
##                          selected_loc_var_pairs = which(loc_var_pairs_partition==i), 
##                          Vecchia_approx_DAG,
##                          y
##                        )
##   )
## }
## 
## field_obs_precisions = get_field_obs_precisions(
##   noise_info, 
##   X_noise_list, 
##   useful_stuff, 
##   loc_var_pairs_partition, 
##   Vecchia_approx_DAG,
##   y
## )





