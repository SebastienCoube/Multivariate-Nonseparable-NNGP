

locs = matrix(runif(1000), ncol = 2)
var_tag = 1 + floor(5 * runif(i)) ; var_tag = match(var_tag, unique(var_tag))
m_whatever_closest = 5
m_same_var = 5
m_other_vars = 2

i = 20
plot(locs[seq(i),], pch = var_tag[seq(i)])
points(locs[i,,drop =F], col = 2, pch = var_tag[i])
points(locs[NNarrays_same_var[[1]][[i]],], col = 3)
points(locs[NNarrays_same_var[[2]][[i]],], col = 4)
points(locs[NNarrays_same_var[[3]][[i]],], col = 5)
points(locs[NNarrays_same_var[[4]][[i]],], col = 6)
points(locs[NNarrays_same_var[[5]][[i]],], col = 7)

# var_tag must be integer and unique(var_tag) must be increasinf

find_ordered_nn_multi <- function(locs, m_whatever_closest, m_same_var, m_other_vars, var_tag, lonlat = FALSE, space_time = FALSE, st_scale = NULL){
  
  # number of locations
  n <- nrow(locs)
  nvar = length(unique(var_tag))
  mult <- 2
  
  # FNN::get.knnx has strange behavior for exact matches
  # so add a small amount of noise to each location
  ee <- min(apply( locs, 2, stats::sd ))
  locs <- locs + matrix( ee*1e-4*stats::rnorm(n*ncol(locs)), n, ncol(locs) ) 
  
  if(lonlat){ # convert lonlattime to xyztime or lonlat to xyz
    lon <- locs[,1]
    lat <- locs[,2]
    lonrad <- lon*2*pi/360
    latrad <- (lat+90)*2*pi/360
    x <- sin(latrad)*cos(lonrad)
    y <- sin(latrad)*sin(lonrad)
    z <- cos(latrad)
    if(space_time){
      time <- locs[,3]
      locs <- cbind(x,y,z,time)
    } else {
      locs <- cbind(x,y,z)
    }
  }
  
  
  if(space_time){ 
    d <- ncol(locs)-1
    if( is.null(st_scale) ){
      randinds <- sample(1:n, min(n,200))
      dvec <- c(fields::rdist( locs[randinds,1:d,drop=FALSE] ))
      dvec <- dvec[ dvec > 0]
      med1 <- mean(dvec)
      dvec <- c(fields::rdist( locs[randinds, d + 1, drop=FALSE] ))
      dvec <- dvec[ dvec > 0]
      med2 <- mean(dvec)
      st_scale <- c(med1,med2)
    }
    locs[ , 1:d] <- locs[ , 1:d]/st_scale[1]
    locs[ , d+1] <- locs[ , d+1]/st_scale[2]
  }
  
  # splitting location according to variable tag
  locs_split_by_var = lapply(unique(var_tag), function(tag)locs[var_tag == tag, ,drop = FALSE])
  
  # number of available parents of each variable
  number_available_parents_by_var = apply(sapply(unique(var_tag), function(tag)var_tag == tag), 2, function(x) (cumsum(x)-1))
  number_available_parents_by_var = pmax(number_available_parents_by_var, 0)

  # getting bounds to end the brute force searches
  maxval_per_variable = apply(number_available_parents_by_var, 2, function(x)match(max(m_other_vars, m_same_var), x)+1)
  required_parents = matrix(m_other_vars, n, nvar);required_parents[cbind(seq(n), var_tag), drop = F] = m_same_var
  available_parents_whatever_closest = number_available_parents_by_var - required_parents; available_parents_whatever_closest = pmax(available_parents_whatever_closest, 0)
  available_parents_whatever_closest  = apply(available_parents_whatever_closest, 1, sum)
  maxval_whatever_closest = max(which(available_parents_whatever_closest <= m_whatever_closest))+1
  

  # to store the nearest neighbor indices
  NNarray <- matrix(NA,n,min(n, m_whatever_closest+ m_same_var+ m_other_vars)+1)
  
  # Start with quotas NNarrays
  # do first by brute force  
  NNarrays_same_var = list()
  if(m_same_var>0 | m_other_vars>0)
  {
    for(i in seq(nvar))
    {
      data_inds = which(var_tag[seq(maxval_per_variable[i])]==i)
      dmat = fields::rdist(locs[data_inds,,drop = F], locs[seq(maxval_per_variable[i]),,drop = F])
      dmat[
      matrix(data_inds, ncol = 1)  %x% matrix(1, 1, maxval_per_variable[i])>=#parent index
      col(dmat)#child index
      ] = Inf # forbidding locations coming after data to be parents
      ordered_dmat = apply(dmat, 2, order)
      NNarrays_same_var[[i]]= lapply(seq(maxval_per_variable[i]), function(j)
      {
        if(var_tag[j]==i) m = m_same_var
        if(var_tag[j]!=i) m = m_other_vars
        return(which(ordered_dmat[,j]<= min(j, m)))
      }
        )
    }
  }
  
  
  
  
  # Whatever Closest NNarray
  # getting first row of whatever closest NNarray
  dmat = fields::rdist(locs[seq(maxval_whatever_closest),,drop = F],)
  dmat[
    row(dmat)>= col(dmat)
  ] = Inf # forbidding locations coming after data to be parents
  for(i in seq(length(NNarrays_same_var))){
    for(j in seq(maxval_whatever_closest))
    {
      dmat[cbind(NNarrays_same_var[[i]][[j]], j)] = Inf
    }
  }
  # forbidding locations already allocated 
  ordered_dmat = apply(dmat, 2, order)

  
  
  
  
  
  
  query_inds <- min( maxval+1, n):n
  data_inds <- 1:n
  
  msearch <- m
  
  while( length(query_inds) > 0 ){
    msearch <- min( max(query_inds), 2*msearch )
    data_inds <- 1:min( max(query_inds), n )
    NN <- FNN::get.knnx( locs[data_inds,,drop=FALSE], locs[query_inds,,drop=FALSE], msearch )$nn.index
    less_than_k <- t(sapply( 1:nrow(NN), function(k) NN[k,] <= query_inds[k]  ))
    sum_less_than_k <- apply(less_than_k,1,sum)
    ind_less_than_k <- which(sum_less_than_k >= m+1)
    
    NN_m <- t(sapply(ind_less_than_k,function(k) NN[k,][less_than_k[k,]][1:(m+1)] ))
    
    NNarray[ query_inds[ind_less_than_k], ] <- NN_m
    
    query_inds <- query_inds[-ind_less_than_k]
    
  }
  
  return(NNarray)
}

