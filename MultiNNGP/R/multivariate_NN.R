###### test data
#####
#####locs_ = matrix(runif(500), ncol = 2)
#####var_tag = 1 + floor(3 * runif(nrow(locs_))) ; var_tag = match(var_tag, unique(var_tag))
#####m_whatever_closest = 10
#####m_same_var = 5
#####m_other_vars = 3
#####
#####
#####tatato = find_ordered_nn_multi(locs_, m_whatever_closest = 0, m_same_var = 5, m_other_vars = 2, var_tag = var_tag)
#####
#####i = 500
#####plot(locs_[seq(i-1),], pch = 16, cex=.5, col = var_tag)
#####points(locs_[i,, drop=F], cex=.5, col = var_tag[i], pch=2)
#####points(locs_[tatato[i,],], col = 4)
#####
##### # test for neighbors from one variable (run inside of function)
##### i = 13
##### plot(locs_[seq(i-1),], pch = 16, cex=.5, col = var_tag)
##### points(locs_[i,, drop = F], cex=1, col = var_tag[i])
##### points(locs_[NNlists_quotas[[i]],], cex=1, col = var_tag[NNlists_quotas[[i]]])
##### points(locs_[NNlist_closest[[i]],], cex=1, col = 3)
##### 
##### 
##### # test for stratified neighbors
##### NNarray = find_ordered_nn_multi(locs_, m_whatever_closest = 10, m_same_var = 5, m_other_vars = 2, var_tag = var_tag, locs_ = locs_)
##### i = 4
##### plot(locs_[seq(i-1),], pch = 16, cex=.5, col = var_tag)
##### points(locs_[i,, drop = F], cex=1, col = var_tag[i])
##### points(locs_[NNlist[[i]],], cex=1, col = 3)
##### 
##### 
##### # wierd test data
##### locs_ = matrix(runif(500), ncol = 2)
##### var_tag = 1 + (locs_[,1]<.5);var_tag = 1 + (var_tag!=var_tag[1])
##### m_whatever_closest = 10
##### m_same_var = 3
##### m_other_vars = 2
##### 
##### # test for stratified neighbors
##### NNlist = find_ordered_nn_multi(locs_, m_whatever_closest = 0, m_same_var = 5, m_other_vars = 2, var_tag = var_tag, locs_ = locs_)
##### i = 200
##### plot(locs_[seq(i-1),], pch = 16, cex=.5, col = var_tag)
##### points(locs_[i,, drop = F], cex=1, col = var_tag[i])
##### points(locs_[NNlist[[i]],], cex=1, col = 3)
 
# var_tag must be integer and unique(var_tag) must be increasing
find_ordered_nn_multi <- function(locs_, m_whatever_closest, m_same_var, m_other_vars, var_tag, lonlat = FALSE){
  var_tag = match(var_tag, unique(var_tag)) # making unique(var_tag) increasing
  if(m_other_vars==0 & m_same_var==0)
  {
    NNarray = GpGp::find_ordered_nn(locs = locs_, m = m_whatever_closest, lonlat = lonlat)
    NNarray = NNarray[,-1]
    NNarray = split(NNarray, row(NNarray))
    NNarray = lapply(NNarray, na.omit)
    return(NNarray)
  }
  # number of locations
  n <- nrow(locs_)
  nvar = length(unique(var_tag))
  
  # FNN::get.knnx has strange behavior for exact matches
  # so add a small amount of noise to each location
  ee <- apply( locs_, 2, stats::sd )
  ee = ee[ee!=0]
  ee = min(ee)
  locs_ <- locs_ + matrix( ee*1e-4*stats::rnorm(n*ncol(locs_)), n, ncol(locs_) ) 
  
  # number of available parents of each variable
  number_available_parents_by_var = apply(sapply(unique(var_tag), function(tag)var_tag == tag), 2, function(x) (cumsum(c(0, x[-length(x)]))))
  number_available_parents_by_var = pmax(number_available_parents_by_var, 0)

  # first search for quotas NNs
  # then search for whatever closest NNs
  # number of parents for each quota
  required_parents_by_var = matrix(m_other_vars, n, nvar);required_parents_by_var[cbind(seq(n), var_tag)] = m_same_var
  obtained_parents_by_var = pmin(required_parents_by_var, number_available_parents_by_var)
   # number of parents for whatever closest
  obtained_parents_whatever_closest = seq(n)-1 - apply(required_parents_by_var, 1, sum)
  obtained_parents_whatever_closest = pmin(obtained_parents_whatever_closest, m_whatever_closest)
  obtained_parents_whatever_closest = pmax(obtained_parents_whatever_closest, 0)
   # number of available parents left once the quotas are filled
  available_parents_whatever_closest = number_available_parents_by_var - required_parents_by_var; available_parents_whatever_closest = pmax(available_parents_whatever_closest, 0)
  available_parents_whatever_closest  = apply(available_parents_whatever_closest, 1, sum)
   # checking when the available closest parents are more than the requirement
  
  # Start with quotas NNlists
  NNlists_quotas = list()
  for(i in seq(nvar))
  {
    NNlists_quotas[[i]] = list()
    #  Tree search
    query_inds = 1:n # start with queries ranging from maxval to n
    data_inds = which(var_tag==i) # start with data being all indices of the variable
    msearch <- max(m_other_vars, m_same_var) # search depht
    while( length(query_inds) > 0 ){
      data_inds <-  which(var_tag[1:max(query_inds)]==i)# restrict data to all indices of the variable coming before last query
      msearch <- min(length(data_inds), 2*msearch)# increase search depht
      NN <- FNN::get.knnx( locs_[data_inds,,drop=FALSE], locs_[query_inds,,drop=FALSE], msearch )$nn.index # get NNs
      NN <- matrix(apply(NN, 2, function(x)data_inds[x]), ncol = ncol(NN))
      # forbidding points of variable i to take themselves as parents
      NN[which(var_tag[query_inds]==i),1] = n+1
      less_than_k <- t(sapply( 1:nrow(NN), function(k) NN[k,] <= query_inds[k]  ))# get only NNs coming before
      sum_less_than_k <- apply(less_than_k,1,sum)
      m = obtained_parents_by_var[query_inds,i]
      ind_less_than_k <- which(sum_less_than_k >= m)
      NN_m <- t(lapply(ind_less_than_k,function(k) NN[k,][less_than_k[k,]][seq(0, m[k])]))
      NNlists_quotas[[i]][ query_inds[ind_less_than_k] ] = NN_m
      query_inds <- query_inds[-ind_less_than_k]
    }
  }
  #merging NNlists quotas
  for(i in seq(length(NNlists_quotas[[1]])))
  {
    NNlists_quotas[[1]][[i]]= unlist(lapply(NNlists_quotas, function(x)x[[i]]))
  }
  NNlists_quotas = NNlists_quotas[[1]]
  if(m_whatever_closest!=0){
    NNarray_whatever_closest = GpGp::find_ordered_nn(locs_, m = m_same_var +(nvar-1)*m_other_vars + m_whatever_closest)
    NNarray_whatever_closest = NNarray_whatever_closest[,-1]
    NNlist_closest = split(NNarray_whatever_closest, row(NNarray_whatever_closest))
  }
  NNlist_closest = mapply(setdiff, NNlist_closest, NNlists_quotas)
  NNlist_closest = lapply(NNlist_closest, function(x)x[seq(m_whatever_closest)])
  NNlist = mapply(c, NNlist_closest, NNlists_quotas)
  NNlist = lapply(NNlist, function(x)c(na.omit(x)))
  return(NNlist)
}

# var_tag must be integer
find_unordered_nn_multi <- function(locs_, m_whatever_closest, m_same_var, m_other_vars, var_tag, lonlat = FALSE){
  var_tag = match(var_tag, unique(var_tag)) # making unique(var_tag) increasing
  # number of locations
  n <- nrow(locs_)
  nvar = length(unique(var_tag))
  
  # FNN::get.knnx has strange behavior for exact matches
  # so add a small amount of noise to each location
  ee <- min(apply( locs_, 2, stats::sd ))
  locs_ <- locs_ + matrix( ee*1e-4*stats::rnorm(n*ncol(locs_)), n, ncol(locs_) ) 
  
  # Start with quotas NNlists
  NNlists_quotas = list()
  for(i in seq(nvar))
  {
    data_inds <-  which(var_tag==i)# restrict data to all indices of the variable coming before last query
    # getting k nearest neighbors
    NNlists_quotas[[i]] = FNN::get.knnx(data = locs_[data_inds,], query = locs_, k = max(m_other_vars, m_same_var)+1)$nn.index
    # traducig k nn into indices of pooled variables instead of indices of
    NNlists_quotas[[i]][]= data_inds[NNlists_quotas[[i]]]
    # traducing neighbor index from 
    # splitting by row and removing superfluous neighbors
    NNlists_quotas[[i]] = split(NNlists_quotas[[i]], row(NNlists_quotas[[i]]))
    NNlists_quotas[[i]] = mapply(
      function(NN_vec, tag)
      {
        if(tag==i)if(m_same_var  ==0){return(NULL)}else {return(NN_vec[2:(m_same_var+1)])}# start by 2 to forbid taking oneself as neighbor
        if(tag!=i)if(m_other_vars==0){return(NULL)}else {return(NN_vec[1:(m_other_vars)])}
      },
      NNlists_quotas[[i]], 
      var_tag, SIMPLIFY = F
      )
  }
  
  #merging NNlists quotas
  NNlists_quotas = lapply(lapply(seq(length(NNlists_quotas[[1]])), function(i)lapply(NNlists_quotas, function(x)x[[i]])), unlist)
  
  # NNlist of whatever nearest neighbor
  # get NNs
  NNlist_closest = FNN::get.knn(data = locs_, k = m_whatever_closest + m_same_var + (nvar-1)*m_other_vars)$nn.index
  # split
  NNlist_closest = split(NNlist_closest, row(NNlist_closest))
  # remove already chosen nn because of quotas
  NNlist_closest = mapply(setdiff, NNlist_closest, NNlists_quotas, SIMPLIFY = F)
  # take m_whatever_closest first neighbors
  NNlist_closest = lapply(NNlist_closest, function(x)x[1:m_whatever_closest])
  
  # merging quotas and closest
  NNlist = mapply(c, NNlist_closest, NNlists_quotas, SIMPLIFY = F)
  
  # adding same index because i conditions on i at previous times
  NNlist = mapply(c, lapply(seq(n), function(x)x), NNlist, SIMPLIFY = F)
  NNlist = lapply(NNlist, FUN = unique)
  return(NNlist)
}

#### # test data
#### 
#### locs_ = matrix(runif(500), ncol = 2)
#### var_tag = 1 + floor(3 * runif(nrow(locs_))) ; var_tag = match(var_tag, unique(var_tag))
#### m_whatever_closest = 10
#### m_same_var = 5
#### m_other_vars = 3
#### 
#### 
#### nnlist_ordered = find_ordered_nn_multi(locs_, m_whatever_closest = 5, m_same_var = 0, m_other_vars = 0, var_tag = var_tag)
#### 
#### nnlist_unordered = find_unordered_nn_multi(locs_, m_whatever_closest = 5, m_same_var = 0, m_other_vars = 0, var_tag = var_tag)
#### 
#### i= 40
#### 
#### # same time neighbors
#### plot(locs_[1:i,], col = var_tag, pch = 16, cex = .5)
#### points(locs_[i,,drop = F], col = var_tag[i], pch = 2)
#### points(locs_[i,,drop = F], col = var_tag[i], pch = 2, cex = 2)
#### points(locs_[nnlist_ordered[[i]],], col = var_tag[nnlist_ordered[[i]]], pch = 1)
#### 
#### 
#### # previous neighbors
#### plot(locs_, col = var_tag, pch = 16, cex = .5)
#### points(locs_[i,,drop = F], col = var_tag[i], pch = 2, cex = 2)
#### points(locs_[nnlist_unordered[[i]],], col = var_tag[nnlist_unordered[[i]]], pch = 1)






make_simple_Vecchia_approx_DAG = 
  function(
    y, 
    locs, 
    m_same_var_same_time = 5, 
    m_other_vars_same_time = 0, 
    m_whatever_closest_same_time = 10, 
    m_same_var_previous_times = 5, 
    m_other_vars_previous_times = 0, 
    m_whatever_closest_previous_times = 10
  ){
    # splitting y and loc in a loc-var format
    var_tag = t(outer(rep(1, dim(y)[1]), seq(dim(y)[2])))
    loc_idx = t(outer(seq(dim(y)[1]), rep(1, dim(y)[2])))
    at_least_one_obs = t(!apply(y, c(1, 2), function(x)all(is.na(x)))) # loc - var pairs with at least one obs
    locs_ = locs[loc_idx[at_least_one_obs],]# expanding
    var_tag_ = var_tag[at_least_one_obs]
    #
    NNarray_same_time = find_ordered_nn_multi(
      locs = locs_, 
      var_tag = var_tag_, 
      m_whatever_closest = m_whatever_closest_same_time,
      m_same_var         = m_same_var_same_time,
      m_other_vars       = m_other_vars_same_time,
      lonlat = F)
    NNarray_pevious_times = NULL
    if(dim(y)[3]>1)NNarray_pevious_times = find_unordered_nn_multi(
      locs = locs_, 
      var_tag = var_tag_, 
      m_whatever_closest = m_whatever_closest_previous_times, 
      m_same_var = m_same_var_previous_times, 
      m_other_vars = m_other_vars_previous_times, 
      lonlat = F)
    Vecchia_approx_DAG = 
      (
        list(
          "DAG" = list(
            "children" = lapply(seq(nrow(locs_)), function(x)x), 
            "parents_same_time" = NNarray_same_time, 
            "parents_previous_times" = NNarray_pevious_times 
          ), 
          field_position = list( # position of sampled w in the loc-var array 
            "location_idx" = loc_idx[at_least_one_obs],
            "var_idx" = var_tag[at_least_one_obs], 
            "loc_match" = split(seq(length(loc_idx[at_least_one_obs])), loc_idx[at_least_one_obs])
          )
        )
      )
    Vecchia_approx_DAG$utils = list(
      lower_tri_idx_DAG = get_lower_tri_idx_DAG(Vecchia_approx_DAG$DAG),
      var_idx = get_var_idx(Vecchia_approx_DAG$DAG, Vecchia_approx_DAG$field_position$var_idx)
    )
    return(Vecchia_approx_DAG)
  }
