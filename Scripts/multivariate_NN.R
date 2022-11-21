

# test data
## locs = matrix(runif(50000), ncol = 2)
## var_tag = 1 + floor(3 * runif(nrow(locs))) ; var_tag = match(var_tag, unique(var_tag))
## m_whatever_closest = 10
## m_same_var = 5
## m_other_vars = 0
## 
## 
## tatato = find_ordered_nn_multi(locs, m_whatever_closest = 0, m_same_var = 5, m_other_vars = 0, var_tag = var_tag)
## 
## i = 1000
## plot(locs[seq(i-1),], pch = 16, cex=.5, col = var_tag)
## points(locs[i,, drop=F], cex=.5, col = var_tag[i], pch=2)
## points(locs[tatato[i,],], col = 4)

#### # test for neighbors from one variable (run inside of function)
#### i = 13
#### plot(locs[seq(i-1),], pch = 16, cex=.5, col = var_tag)
#### points(locs[i,, drop = F], cex=1, col = var_tag[i])
#### points(locs[NNlists_quotas[[i]],], cex=1, col = var_tag[NNlists_quotas[[i]]])
#### points(locs[NNlist_closest[[i]],], cex=1, col = 3)
#### 
#### 
#### # test for stratified neighbors
#### NNarray = find_ordered_nn_multi(locs, m_whatever_closest = 10, m_same_var = 5, m_other_vars = 2, var_tag = var_tag, locs = locs)
#### i = 4
#### plot(locs[seq(i-1),], pch = 16, cex=.5, col = var_tag)
#### points(locs[i,, drop = F], cex=1, col = var_tag[i])
#### points(locs[NNlist[[i]],], cex=1, col = 3)
#### 
#### 
#### # wierd test data
#### locs = matrix(runif(500), ncol = 2)
#### var_tag = 1 + (locs[,1]<.5);var_tag = 1 + (var_tag!=var_tag[1])
#### m_whatever_closest = 10
#### m_same_var = 3
#### m_other_vars = 2
#### 
#### # test for stratified neighbors
#### NNlist = find_ordered_nn_multi(locs, m_whatever_closest = 0, m_same_var = 5, m_other_vars = 2, var_tag = var_tag, locs = locs)
#### i = 200
#### plot(locs[seq(i-1),], pch = 16, cex=.5, col = var_tag)
#### points(locs[i,, drop = F], cex=1, col = var_tag[i])
#### points(locs[NNlist[[i]],], cex=1, col = 3)
#### 
# var_tag must be integer and unique(var_tag) must be increasing
find_ordered_nn_multi <- function(locs, m_whatever_closest, m_same_var, m_other_vars, var_tag, lonlat = FALSE){
  if(m_other_vars==0 & m_same_var==0)return(GpGp::find_ordered_nn(locs = locs, m = m_whatever_closest, lonlat = lonlat))
  # number of locations
  n <- nrow(locs)
  nvar = length(unique(var_tag))
  mult <- 2
  
  # FNN::get.knnx has strange behavior for exact matches
  # so add a small amount of noise to each location
  ee <- min(apply( locs, 2, stats::sd ))
  locs <- locs + matrix( ee*1e-4*stats::rnorm(n*ncol(locs)), n, ncol(locs) ) 
  
  
  # splitting location according to variable tag
  locs_split_by_var = lapply(unique(var_tag), function(tag)locs[var_tag == tag, ,drop = FALSE])
  
  # number of available parents of each variable
  number_available_parents_by_var = apply(sapply(unique(var_tag), function(tag)var_tag == tag), 2, function(x) (cumsum(x)-1))
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
      NN <- FNN::get.knnx( locs[data_inds,,drop=FALSE], locs[query_inds,,drop=FALSE], msearch )$nn.index # get NNs
      NN <- apply(NN, 2, function(x)data_inds[x])
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
  NNlist_closest = lapply(seq(n), function(i)NULL)
  if(m_whatever_closest!=0){
    # Whatever closest NNlist
    NNlist_closest = list()
    query_inds = 1:n # start with queries ranging from 1 to n
    data_inds = 1:n # start with data being all indices of the variable
    msearch <- max(m_whatever_closest) # search depht
     while( length(query_inds) > 0 ){
      data_inds <-  1:max(query_inds)# restrict data to all indices of the variable coming before last query
      msearch <- min(length(data_inds), 2*msearch)# increase search depht
      NN <- FNN::get.knnx( locs[data_inds,,drop=FALSE], locs[query_inds,,drop=FALSE], msearch )$nn.index # get NNs
      # forbidding parents already chosen for quotas
      already_taken_indices = mapply(function(x, y) which(x%in%y), split(NN, row(NN)), NNlists_quotas[query_inds], SIMPLIFY = F)
      NN[cbind(unlist(lapply(seq(length(already_taken_indices)), function(i)rep(i, length(already_taken_indices[[i]])))), unlist(already_taken_indices))] = n+1
      NN[,1] = n+1
      less_than_k <- t(sapply( 1:nrow(NN), function(k) NN[k,] <= query_inds[k]  ))# get only NNs coming before
      sum_less_than_k <- apply(less_than_k,1,sum)
      m = obtained_parents_whatever_closest[query_inds]
      ind_less_than_k <- which(sum_less_than_k >= m )
      NN_m <- t(lapply(ind_less_than_k,function(k) NN[k,][less_than_k[k,]][1:(m[k])]))
      NNlist_closest[ query_inds[ind_less_than_k] ] <- NN_m
      query_inds <- query_inds[-ind_less_than_k]
    }
  }
  NNlist = mapply(c, NNlist_closest, NNlists_quotas)
  NNlist = lapply(NNlist, function(x)c(na.omit(x)))

  NNarray = matrix(NA, length(NNlist), max(sapply(NNlist, length))+1)
  ##NNarray[cbind(
  ##  rep(seq(length(NNlist)), sapply(NNlist, length)), 
  ##  unlist(lapply(NNlist, function(x)seq(length(x))))[-c(1, 2)]+1
  ##  )] = unlist(NNlist)
  NNarray[,-1][
    cbind(
      rep(seq(length(NNlist)), sapply(NNlist, length)),
      unlist(lapply(NNlist, function(x)
      {
        if(length(x)==0)return(c())
        if(length(x)!=0)return(seq(length(x)))
      })))] = unlist(NNlist)
  NNarray[,1] = seq(n)
  return(NNarray)
}
