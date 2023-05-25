# To do: improve class into tiles and get basis functions, 
# can be much more efficient 



get_blocks_n_bases_with_witnesses = 
  function(mcmc_nngp_list, 
           space_cluster_size_target = 100){
  ###########################
  # SPATIAL BASIS FUNCTIONS #
  ###########################
  basis_functions = c(
    # groups of spatial basis functions
    unlist(lapply(
      get_grids(
        points = mcmc_nngp_list$useful_stuff$locs_repeated, 
        cluster_size_target = space_cluster_size_target), 
      function(x)get_basis_functions(
        points = mcmc_nngp_list$useful_stuff$locs_repeated, 
        tile_info = x, 
        cluster_size_target = space_cluster_size_target, 
        Vecchia_approx_DAG = mcmc_nngp_list$Vecchia_approx_DAG, 
        useful_stuff = mcmc_nngp_list$useful_stuff)),
      recursive = F)
  )
  # partitioning the space into clusters whose size is lower than cluster_size_target
  recursive_k_means_clust = recursive_k_means(locs_ = mcmc_nngp_list$useful_stuff$locs_repeated, space_cluster_size_target)
  # removing indicator basis functions to guarantee reversibility
  # by construction the basis functions are independent
  # removing the lowest index in order to reduce the number of colors
  # as many removed indices as basis functions
  removed_indices = rep(0, sum(sapply(basis_functions, ncol)))
  candidate_loc_idx = seq(nrow(basis_functions[[1]]))
  count = 1
  # chosing witnesses of smallest order in Vecchia ordering
  for(i_basis_function in seq(length(basis_functions)))
  {
    for(j in seq(ncol(basis_functions[[i_basis_function]])))
    {
      # choosing smallest index available among where basis function is != 0
      removed_indices[count]= min(
        candidate_loc_idx[which(basis_functions[[i_basis_function]][,j] != 0)]
      )
      # forbidding index to be chosen again
      candidate_loc_idx[which(basis_functions[[i_basis_function]][,j] != 0)
      ][which.min(
        candidate_loc_idx[which(basis_functions[[i_basis_function]][,j] != 0)]
      )] = Inf
      count = count+1
    }
  }
  # test : does any basis function not have a witness ?
  ## any(
  ##   (
  ##     Matrix::sparseVector(x = 1, i = removed_indices, length = nrow(basis_functions[[1]])) %*%
  ##       (do.call(cbind, basis_functions)!=0)
  ##   )==0)
  recursive_k_means_clust$clust[removed_indices] = NA
  # adding blocks to bases
  basis_functions = 
    c(
      basis_functions, 
      get_indicator_basis_functions (
        recursive_k_means_clust = recursive_k_means_clust, 
        useful_stuff = mcmc_nngp_list$useful_stuff)
    )
  spatial_coloring = 
    color_basis_functions(
      basis_functions = basis_functions, 
      precision_blocks = chain$stuff$precision_blocks
    )
  ############
  # COLORING #
  ############
  return(list(
    basis_functions = basis_functions, 
    spatial_coloring = spatial_coloring
    ))
}




# split vectors of time with some randomness of frontier
split_vector <- function(vec, target_size = 100) {
  n <- length(vec)
  k <- max(floor(n/target_size), 1)  # determine number of subvectors needed
  m = ceiling(n/k)
  subvecs <- vector("list", k)
  start <- 1
  for (i in 1:k) {
    end <- min(i*m + rbinom(1, 1, .5)*min(rpois(1,10), 15), n)
    subvecs[[i]] <- vec[start:end]
    start <- end + 1
  }
  return(
    cbind(
      sapply(subvecs, function(x)x[1]),
      sapply(subvecs, function(x)x[length(x)])
    )
  )
}


get_time_split = function(time_depth, n_time_periods, time_target_size)
{
  time_split = split_vector(
    seq(
      time_depth, 
      n_time_periods - time_depth +1
    ), 
    target_size = time_target_size
  )
  list(time_split = time_split, time_coloring = rep(seq(nrow(time_split))%%2))
}

recursive_k_means = function(locs_, cluster_size_target)
{
  clust = rep(1, nrow(locs_))
  while(any(table(clust)>cluster_size_target))
  {
    for(i in 6)
    {
      for(j in which(table(clust)>(i-1)*cluster_size_target))
      {
        clust[clust == j] = 
          paste(clust[clust == j], kmeans(locs_[clust == j,], centers = 3+ceiling(3*runif(1)))$cl, sep = "_")
      }
      clust = match(clust, unique(clust))
    }
    for(i in seq(3, 2))
    {
      for(j in which(table(clust)>(i-1)*cluster_size_target))
      {
        clust[clust == j] = 
          paste(clust[clust == j], kmeans(locs_[clust == j,], centers = i)$cl, sep = "_")
      }
      clust = match(clust, unique(clust))
    }
  }
  res = list()
  res$centers = t(sapply(unique(clust), function(i)apply(locs_[clust==i,,drop =F], 2, mean)))
  res$clust = clust
  res
}

get_indicator_basis_functions = function(recursive_k_means_clust, useful_stuff)
{
  lapply(as.vector(na.omit(unique(recursive_k_means_clust$clust))), function(i)
    Matrix::sparseMatrix(
      i = which(recursive_k_means_clust$clust==i),
      j = seq(sum(na.omit(recursive_k_means_clust$clust==i))),
      x=1, dims = c(useful_stuff$n_field, sum(na.omit(recursive_k_means_clust$clust==i)))
    )
  )
}
## get_indicator_basis_functions_REMOVE_HIGHEST = function(recursive_k_means_clust, useful_stuff)
## {
##   lapply(unique(recursive_k_means_clust$clust), function(i)
##   {
##     idx = which(recursive_k_means_clust$clust==i)
##     idx = idx[-seq(5)]
##     Matrix::sparseMatrix(
##       i = idx,
##       j = seq(length(idx)),
##       x=1, dims = c(useful_stuff$n_field, length(idx))
##     )
##   }
##   )
## }

color_basis_functions_time_slicing = function(basis_functions, precision_blocks, time_markov_mat)
{
  basis_functions_nonzero_sites = 
    do.call(
      cbind, 
      lapply(basis_functions, 
             function(x)x%*%Matrix::sparseVector(x = rep(1, ncol(x)), i = seq(ncol(x)), length = (ncol(x))))
    )
  basis_functions_markov_mat = Matrix::t(basis_functions_nonzero_sites) %*% (precision_blocks[[1]]+precision_blocks[[2]]) %*% basis_functions_nonzero_sites
  basis_functions_markov_mat@x = rep(1, length(basis_functions_markov_mat@x))
  markov_mat = time_markov_mat %x% basis_functions_markov_mat
  naive_greedy_coloring(markov_mat)
}

color_basis_functions = function(basis_functions, precision_blocks)
{
    basis_functions_nonzero_sites = 
    do.call(
      cbind, 
      lapply(basis_functions, 
             function(x)x%*%Matrix::sparseVector(x = rep(1, ncol(x)), i = seq(ncol(x)), length = (ncol(x))))
    )
  basis_functions_markov_mat = Matrix::t(basis_functions_nonzero_sites) %*% (precision_blocks[[1]]+precision_blocks[[2]]) %*% basis_functions_nonzero_sites
  basis_functions_markov_mat@x = rep(1, length(basis_functions_markov_mat@x))
  DSATUR(basis_functions_markov_mat)
}


# Define the function to class the locs into square tiles
class_points_into_tiles <- function(locs, tile_size) {
  # Compute the minimum and maximum x and y coordinates of the locs
  min_x <- min(locs[,1])
  max_x <- max(locs[,1])
  min_y <- min(locs[,2])
  max_y <- max(locs[,2])
  
  # Compute the number of tiles in each dimension
  n_x_tiles <- ceiling((max_x - min_x) / tile_size)+1
  n_y_tiles <- ceiling((max_y - min_y) / tile_size)+1
  
  # Classify each point into a tile
  #tile_indices <- matrix(0, nrow = nrow(locs), ncol = 2)
  #for (i in 1:nrow(locs)) {
  #  x_tile <- ceiling((locs[i,1] - min_x +1e-8) / tile_size)
  #  y_tile <- ceiling((locs[i,2] - min_y +1e-8) / tile_size)
  #  tile_indices[i,] <- c(x_tile, y_tile)
  #}
  
  tile_indices = cbind(
    as.numeric(cut(locs[,1]+abs(rnorm(nrow(locs), 0, 1e-10)), breaks = unique(c(seq(min_x-.00001, max_x + .00001, tile_size), max_x)))), 
    as.numeric(cut(locs[,2]+abs(rnorm(nrow(locs), 0, 1e-10)), breaks = unique(c(seq(min_y-.00001, max_y + .00001, tile_size), max_x))))
  )
  
  tile_centers <- array(0, dim = c(n_x_tiles, n_y_tiles, 2))
  #x_noise = tile_size* runif(1)
  #y_noise = tile_size* runif(1)
  for (i in 1:n_x_tiles) {
    for (j in 1:n_y_tiles) {
      x_center <- min_x + (i - 0.5) * tile_size# - x_noise
      y_center <- min_y + (j - 0.5) * tile_size# - y_noise
      tile_centers[i,j,] <- c(x_center, y_center)
    }
  }
  # Return the tile indices and centers
  return(list(tile_indices = tile_indices, tile_centers = tile_centers, tile_size = tile_size))
}

## Example usage
#set.seed(123)
#points <- matrix(runif(50000), ncol = 2)
#tile_size <- 0.5
#
## Classify the points into tiles and compute the centers
#tile_info <- class_points_into_tiles(points, tile_size)


get_grids = function(points, cluster_size_target)
{
  d_max = max(apply(points, 2, max)- apply(points, 2, min))
  tile_size = d_max/11.9
  res = list()
  tile_info = class_points_into_tiles(points, tile_size)
  res[[1]] = tile_info
  tatato = 1
  while(any(table(tile_info$tile_indices[,1], tile_info$tile_indices[,2])>ceiling(sqrt(cluster_size_target)))){
    tile_size = tile_size/ceiling(sqrt(cluster_size_target))
    tile_info = class_points_into_tiles(points, tile_size)
    tatato = tatato +1
    res[[tatato]] = tile_info
  }
  res[length(res)]=NULL
  #plot(points, col = interaction(tile_info$tile_indices[,1], tile_info$tile_indices[,2]))
  #table(interaction(tile_info$tile_indices[,1], tile_info$tile_indices[,2]))
  return(res)
}

#points = locs_
#cluster_size_target = 30
#grids = get_grids(points = points, cluster_size_target = cluster_size_target)
#table(grids[[1]]$tile_indices[,1], grids[[1]]$tile_indices[,2])
#table(grids[[2]]$tile_indices[,1], grids[[2]]$tile_indices[,2])
#table(grids[[3]]$tile_indices[,1], grids[[3]]$tile_indices[,2])
#table(grids[[9]]$tile_indices[,1], grids[[9]]$tile_indices[,2])
#
#plot(points, col = interaction(grids[[3]]$tile_indices[,1], grids[[3]]$tile_indices[,2]))
#
#tile_info = grids[[9]]


get_basis_functions = function(points, tile_info, cluster_size_target, Vecchia_approx_DAG, useful_stuff)
{
  size_table = table(tile_info$tile_indices[,1], tile_info$tile_indices[,2])>ceiling(sqrt(cluster_size_target))
  active_centers = list()
  i_idx = list()
  j_idx = list()
  x = list()
  # looping on the grid
  for(center_i in seq(nrow(tile_info$tile_centers)))
  {
    for(center_j in seq(ncol(tile_info$tile_centers)))
    {
      # selecting tiles of the grid with more than cluster_size_target points in it
      if(!is.na(size_table[match(center_i, row.names(size_table)),match(center_j, colnames(size_table))])){
        if(size_table[match(center_i, row.names(size_table)),match(center_j, colnames(size_table))])
        {
          # appending current center to active centers
          active_centers = c(active_centers, list(tile_info$tile_centers[center_i, center_j,]))
          # selecting locs in adjacent tiles
          selected_locs = 
            which(
              (center_i-3 < tile_info$tile_indices[,1]) &
              (center_i+3 > tile_info$tile_indices[,1]) &
              (center_j-3 < tile_info$tile_indices[,2]) &
              (center_j+3 > tile_info$tile_indices[,2]) 
            )
          # indices of spatial basis function
          i_idx = c(i_idx, list(selected_locs))
          x = c(x, list(c(exp(-(fields::rdist(
            x1 = matrix(tile_info$tile_centers[center_i, center_j,], ncol = 2), 
            x2 = points[selected_locs,]
            )^2/tile_info$tile_size^2)))))
        }
      }
    }
  }
  
  # making basis functions
  basis_functions = mapply(
    function(ii, xx)Matrix::sparseMatrix(
      i = ii, 
      j = Vecchia_approx_DAG$field_position$var_idx[ii], 
      x = xx, dims = c(useful_stuff$n_field, useful_stuff$n_var_y)
    ), 
    ii = i_idx, 
    xx = x
  )
  # removing all  zero columns
  basis_functions = lapply(basis_functions, function(x)x[,which(as.vector(rep(1, nrow(x)) %*% x) !=0)])
  
  
  center_split = recursive_k_means(
    locs_ = do.call(rbind, active_centers)[rep(seq(length(basis_functions)), sapply(basis_functions, ncol)),], 
    cluster_size_target = cluster_size_target)
  basis_functions_cl = center_split$clust [!duplicated(rep(seq(length(basis_functions)), sapply(basis_functions, ncol)))]
  
  lapply(unique(basis_functions_cl), function(cl)
    do.call(cbind, basis_functions[basis_functions_cl==cl])
    )
}


#spatial_basis_functions = do.call(c, lapply(grids, function(g)get_basis_functions(points = points, tile_info = g, cluster_size_target = cluster_size_target, Vecchia_approx_DAG = Vecchia_approx_DAG, useful_stuff = useful_stuff)))
#plot(locs_[,1], spatial_basis_functions[[73]][,1])
#plot(locs_[,1], spatial_basis_functions[[20]][,8])
#plot(locs_[,1], spatial_basis_functions[[1]][,8])
#plot(locs_[,1], spatial_basis_functions[[2]][,8])
#plot(locs_[,1], spatial_basis_functions[[4]][,8])
#plot(locs_[,1], spatial_basis_functions[[5]][,8])
#plot(locs_[,1], spatial_basis_functions[[9]][,8])
#plot(locs_[,1], spatial_basis_functions[[9]][,14])


