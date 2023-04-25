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
  lapply(unique(recursive_k_means_clust$clust), function(i)
    Matrix::sparseMatrix(
      i = which(recursive_k_means_clust$clust==i),
      j = seq(sum(recursive_k_means_clust$clust==i)),
      x=1, dims = c(useful_stuff$n_field, sum(recursive_k_means_clust$clust==i))
    )
  )
}

color_basis_functions = function(basis_functions)
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
  tile_indices <- matrix(0, nrow = nrow(locs), ncol = 2)
  for (i in 1:nrow(locs)) {
    x_tile <- ceiling((locs[i,1] - min_x +1e-8) / tile_size)
    y_tile <- ceiling((locs[i,2] - min_y +1e-8) / tile_size)
    tile_indices[i,] <- c(x_tile, y_tile)
  }
  
  tile_centers <- array(0, dim = c(n_x_tiles, n_y_tiles, 2))
  x_noise = tile_size* runif(1)
  y_noise = tile_size* runif(1)
  for (i in 1:n_x_tiles) {
    for (j in 1:n_y_tiles) {
      x_center <- min_x + (i - 0.5) * tile_size - x_noise
      y_center <- min_y + (j - 0.5) * tile_size - y_noise
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


