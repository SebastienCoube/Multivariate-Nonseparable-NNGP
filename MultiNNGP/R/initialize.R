### source("MultiNNGP/R/GMA.R")
### source("MultiNNGP/R/multivariate_NN.R")
### source("MultiNNGP/R/vecchia.R")
 my_image = function(m)image(t(m)[,nrow(m):1])


## # generating data 
## n_loc = 500
## n_time = 200
## n_var = 5
## 
## # full data
## y = array(dim = c(n_loc, n_var, n_time))
## dimnames(y) = list(
##   paste("loc", seq(n_loc), sep = "_"), 
##   paste("var", seq(n_var), sep = "_"), 
##   paste("time", seq(n_time), sep = "_")
##   )
## y[] = rnorm(length(y))
## 
## X = array(dim = c(n_loc, 10, n_time))
## dimnames(X) = list(
##   paste("loc", seq(n_loc), sep = "_"), 
##   paste("covariate", seq(dim(X)[2]), sep = "_"), 
##   paste("time", seq(n_time), sep = "_")
##   )
## X[] = rnorm(length(X))
## X[,1,] = 1
## X_noise = X
## X_scale = X
## 
## locs = matrix(runif(2*n_loc), n_loc)

# computes cross-product of array along the first and third dimension (typically gives a crossprod per variable)
array_crossprod = function(array_1, array_2 = NULL)
{
  if(is.null(array_2))
  {
    res =  matrix(0, dim(array_1)[2], dim(array_1)[2])
    for(i in seq(ncol(res)))
    {
      for(j in seq(1, i))
      {
        res[i, j] = sum(na.omit(c(array_1[,i,]*array_1[,j,])))
        res[j, i] = res[i, j]
      }
    }
  }
  if(!is.null(array_2))
  {
    res =  matrix(0, dim(array_1)[2], dim(array_2)[2])
    for(i in seq(nrow(res)))
    {
      for(j in seq(ncol(res)))
      {
        res[i, j] = sum(na.omit(c(array_1[,i,]*array_2[,j,])))
      }
    }
  }
  res
}

# computes product of array and matrix following the second dimension of the array and the first dimension of the matrix
# array-matrix equivalent of matrix-vector X %*% beta
# a of size j,k,l,   m of size k,p, a*m of size j,p,l
array_matrix_mult = function(a, m)
{
  res = array(0, dim = c(dim(a)[1], dim(m)[2], dim(a)[3]))
  for(i in seq(nrow(m)))for(j in seq(ncol(m)))res[,j,] = res[,j,] + a[,i,]*m[i,j]
  res
}


#  symmetric matrix
get_symmetric_mat = function(mat_coordinates)
{
  M = matrix(0,nrow_square_matrix(length(mat_coordinates)),nrow_square_matrix(length(mat_coordinates)))
  M[lower.tri(M, diag = T)] = mat_coordinates
  M = t(M)
  M[lower.tri(M, diag = T)] = mat_coordinates
  M
}
nrow_square_matrix = function(mat_coord_length)(-1+ sqrt(1 + 8*mat_coord_length))/2


# get exponential of symmetric matrix and its derivatives 
get_expmat_and_derivatives = function(mat_coordinates)
{ 
 M_0 = get_symmetric_mat(mat_coordinates)
 eig = eigen(M_0)
 # using directional derivative formula
 G = outer(eig$val, eig$val, function(x, y)(exp(x)-exp(y))/(x-y))
 diag(G)=exp(eig$val)
 G[is.nan(G)]=0
 J = get_lower_tri_idx(ncol(M_0), T)
 list(
   "expmat" = t(eig$vectors)%*%(exp(eig$values)*eig$vectors),
   "expmat_derivative" = lapply(seq(length(mat_coordinates)), function(i)
   {
     #coord = rep(0, length(mat_coordinates)); coord[i]=1; V = get_symmetric_mat(coord) 
     #eig$vectors %*% ((t(eig$vec) %*% V %*% eig$vectors)*G)%*% t(eig$vectors)
     cross  = tcrossprod(eig$vectors[J[i,1],], eig$vectors[J[i,2],])
     #cross+t(cross)-t(eig$vec) %*% V %*% eig$vectors
     eig$vectors %*% ((cross+t(cross))*G)%*% t(eig$vectors)
   })
 )
}

get_tau_precision = function(expmat, y)solve(expmat[!is.na(y), !is.na(y)])

# get the noise precision matrix and its derivatives following NA pattern of y 
get_tau_info = function(mat_coordinates, y)
{
  if(all(is.na(y)))return(NULL)
  if(all(is.na(mat_coordinates)))return(NULL)
  res = list()
  expmat_and_derivatives = get_expmat_and_derivatives(mat_coordinates = mat_coordinates)
  res$tau_precision = get_tau_precision(expmat_and_derivatives$expmat, y)
  res$tau_precision_deriv = lapply(expmat_and_derivatives$expmat_derivative, function(x) res$tau_precision %*% x[!is.na(y), !is.na(y)] %*% res$tau_precision)
  return(res)
}


### get all noise precision matrices and their derivatives following NA pattern of y 
##get_noise_info = function(X_noise_list, noise_beta, y, y_NA_possibilities)
##{
##  res = list()
##  mat_coordinates_array = array_matrix_mult(X_noise_list$X_killed_NA, noise_beta)
##  # nonstationary case
##  if(!X_noise_list$is_only_intercept) expmat_and_derivatives = parallel::mcmapply(
##    FUN = get_tau_info, 
##    mat_coordinates = lapply(X_noise_list$no_NA_idx_split, function(x)mat_coordinates_array[x[1],,x[2]]),
##    y = lapply(X_noise_list$no_NA_idx_split, function(x)y[x[1],,x[2]]), 
##    SIMPLIFY = F, mc.cores = parallel::detectCores()-1
##  )
##  # stationary case
##  if(X_noise_list$is_only_intercept)
##  {
##    tau_info_cases = lapply(
##      split(y_NA_possibilities, row(y_NA_possibilities)), 
##      function(x)get_tau_info(mat_coordinates = rep(1, nrow(noise_beta))%*%noise_beta, y = x)
##    )
##    expmat_and_derivatives = lapply(
##      X_noise_list$y_NA_pattern, 
##      function(x)tau_info_cases[[x]])
##  }
##    
##  return(expmat_and_derivatives)
##}


# get all noise precision matrices and their derivatives following NA pattern of y 
get_noise_info = function(X_noise_list, noise_beta, y_NA_possibilities_match, y_NA_possibilities, y)
{
  mat_coordinates_array = array_matrix_mult(X_noise_list$X, noise_beta)
  # nonstationary case
  if(!X_noise_list$is_only_intercept) expmat_and_derivatives = parallel::mclapply(
    seq(dim(y)[1]), function(loc_index)
    {
      lapply(seq(dim(y)[3]), function(time_index)
      {
        return(get_tau_info(
          mat_coordinates = X_noise_list$X[loc_index,,time_index]%*%noise_beta,
          y = y[loc_index,,time_index]
          ))
      })
    }, mc.cores = parallel::detectCores()-1)
  # stationary case
  if(X_noise_list$is_only_intercept)
  {
    tau_info_cases = lapply(
      split(y_NA_possibilities, row(y_NA_possibilities)), 
      function(x)
      {
        x[x==F]=NA
        get_tau_info(mat_coordinates = rep(1, nrow(noise_beta))%*%noise_beta, y = x)
      }
    )
    expmat_and_derivatives = lapply(seq(dim(y)[1]), function(loc_index)
   {
     lapply(seq(dim(y)[3]), function(time_index)
     {
       tau_info_cases[[y_NA_possibilities_match[loc_index, time_index]]]
     })
   })
  }
  return(expmat_and_derivatives)
}


#mat_coordinates = rnorm(15) 
#expmat_and_derivatives = get_expmat_and_derivatives(mat_coordinates)
#y = rnorm(5); y[c(3, 5)] = NA
#
#-.5 * log(det(expmat_and_derivatives$expmat[c(1, 2, 4), c(1, 2, 4)]))-
#  .5 * na.omit(y)%*%solve(expmat_and_derivatives$expmat[c(1, 2, 4), c(1, 2, 4)])%*%na.omit(y)


process_covariates = function(X, y)
{
  res = list()
  res$X = X # X itself
  res$X_killed_NA = X; res$X_killed_NA[is.na(res$X_killed_NA)] = 0 # X where NAs replaced by 0
  res$chol_X = chol(crossprod(apply(res$X_killed_NA, 2, c))) # chol of covariance along dim 2
  res$no_NA_pattern =  !apply(X, c(1, 3), anyNA) # complete observations along dim 2
  #res$no_NA_idx = cbind(row(res$no_NA_pattern)[res$no_NA_pattern], col(res$no_NA_pattern)[res$no_NA_pattern]) # dim 1 and 3 of complete observations 
  #res$no_NA_idx_split = split(res$no_NA_idx, row(res$no_NA_idx))  # dim 1 and 3 of complete observations, split into list
  res$is_only_intercept = (dim(X)[2]==1)&all(na.omit(c(X))==1) # bool to check if X is only an intercept
  
  # used in intercept only. tells which y are observed at the indices of res$no_NA_idx_split
  y_NA_possibilities = as.matrix(expand.grid(lapply(seq(dim(y)[2]), function(i)c(1,NA))))
  y_NA_possibilities = y_NA_possibilities[-nrow(y_NA_possibilities),]
  y_NA_possibilities = split(y_NA_possibilities, row(y_NA_possibilities))
  res$y_NA_pattern = match(
    lapply(res$no_NA_idx_split, function(x)is.na(unname(y[x[1],,x[2]]))[seq(dim(y)[2])]),
    lapply(y_NA_possibilities, is.na)
  )
  res
}


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
   return(
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
   }


### Vecchia_approx_DAG = make_simple_Vecchia_approx_DAG(
###   y = y, locs = locs
### )


# computes vecchia approx from parameters
vecchia_block_approx = function(
    Vecchia_approx_DAG, locs, lower_tri_idx_DAG, var_idx, time_depth, #does not depend on params
    rho_vec, a, b, cc, delta, lambda, r, 
    A_vec, nu_vec, a2_vec
){
  if(is.null(a))a=.5
  if(is.null(b))b=.5
  if(is.null(cc))cc=.5
  if(is.null(delta))delta=.5
  if(is.null(lambda))lambda=.5
  if(is.null(r))r=.5
  if(is.null(A_vec))A_vec = rep(.5, length(nu_vec))
  multiplier = get_multiplier(
    a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
    r = r, A_vec = A_vec, nu_vec = nu_vec, a2_vec = a2_vec, 
    u = seq(0, time_depth-1)
  )
  effective_range = get_effective_range(
    a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
    r = r, A_vec = A_vec, a2_vec = a2_vec, u = seq(0, time_depth-1)
  )
  get_vecchia_blocks(
    DAG = Vecchia_approx_DAG$DAG, 
    coeffs = get_linv_coeffs(
      DAG = Vecchia_approx_DAG$DAG, 
      locs_ = locs[Vecchia_approx_DAG$field_position$location_idx,], 
      lower_tri_idx_DAG = lower_tri_idx_DAG, 
      var_idx = var_idx, 
      var_tag = Vecchia_approx_DAG$field_position$var_idx, 
      multiplier = multiplier, 
      effective_range = effective_range, 
      rho_vec = rho_vec, nu_vec = nu_vec), 
    time_depth = time_depth) 
}


### time_depth = NULL
### n_chains = 2

#' @param y a 3-D array 
#' whose first dim corresponds to spatial locations, 
#' second dim corresponds to variable index, 
#' third dim corresponds to the time. 
#' NAs are allowed, and it is encouraged to put a buffer of NAs in the beginning of the time dimension
#' @param X explains the interest variables. a 3-D array  
#' whose first dim corresponds to spatial locations (just like in "y"), 
#' second dim corresponds to variable index (not the same variable index as in "y" !), 
#' third dim corresponds to the time (just like in "y"). 
#' Each observation made in the dimensions 1 and 3 of y (in English : at given location, given time, at least one variable was observed)
#' must be matched by a complete observation in X_per_obs (in English : at given location, given time, all covariates in X were observed).
#' Whenever possible, e.g. in the case of an intercept or basis functions, 
#' put values of X for all times and all places, this will improve greatly the behavior of the associated regression coefficients.
#' @param X_noise explains the noise variance. 
#' Same conditions as X
#' @param X_scale explains the latent field variance. 
#' Must be no NA, even where there are NAs in y. 
#' That is because the latent field is sampled even where the field is not observed
multivariate_NNGP_initialize = function(
    y, locs, 
    X, 
    X_noise,
    X_scale, 
    time_depth = NULL,
    Vecchia_approx_DAG, 
    n_chains = 2
  )
{
  #################
  # Sanity checks #
  #################
  # checking format
  if(!is.array(y))stop("y should be an array")
  if(!is.array(X))stop("X should be an array")
  if(!is.array(X_scale))stop("X_scale should be an array")
  if(!is.array(X_noise))stop("X_noise should be an array")
  if(!is.matrix(locs))stop("locs should be a matrix")
  # checking dimension
  if(dim(y)[1]!= nrow(locs)) stop("the first dimension of y should be equal to the number of spatial locations")
  if(any(dim(y)[c(1,3)]!=dim(X)[c(1,3)]))stop("dimensions 1 (spatial location) and 3 (time)  of y and X should be the same")
  if(any(dim(y)[c(1,3)]!=dim(X_scale)[c(1,3)]))stop("dimensions 1 (spatial location) and 3 (time) of y and X_scale should be the same")
  if(any(dim(y)[c(1,3)]!=dim(X_noise)[c(1,3)]))stop("dimensions 1 (spatial location) and 3 (time) of y and X_noise should be the same")
  # checking X
  if(any((apply(y, c(1, 3), anyNA)==F)&(apply(X, c(1, 3), anyNA)==T)))stop("X should have no NAs where y has no NAs (but it is possible, and even recommended, for X to have observations even where y has NAs)")
  if(any((apply(y, c(1, 3), anyNA)==F)&(apply(X_noise, c(1, 3), anyNA)==T)))stop("X_noise should have no NAs where y has no NAs")
  if(any(is.na(X_scale)))stop("X_scale can have no NA")
  # checking Vecchia DAG
  if(is.null(Vecchia_approx_DAG$DAG$parents_previous_times) & (dim(y)[3]>1))stop("Vecchia_approx_DAG has been made with option `time_depth = 1` while y has more than one time period")
  ###################
  # some processing #
  ###################
  if(is.null(time_depth)){
    if(dim(y)[3]==1)time_depth = 1
    if(dim(y)[3]> 1)time_depth = 5
  }
  # removing unobserved stuff in Y
  all_na_idx = apply(y, c(1), function(x)all(is.na(x)))
  #########################
  # Processing covariates #
  #########################
  # adding buffer
  if(dim(y)[3]>1){
    X =       abind::abind(array(123456 , c(dim(X)      [c(1,2)], 5*time_depth)), X,       along = 3)
    X_scale = abind::abind(array(123456 , c(dim(X_scale)[c(1,2)], 5*time_depth)), X_scale, along = 3)
    X_noise = abind::abind(array(NA,      c(dim(X_noise)[c(1,2)], 5*time_depth)), X_noise, along = 3)
    y =       abind::abind(array(NA ,     c(dim(y)      [c(1,2)], 5*time_depth)), y,       along = 3)
  }
  covariates = 
    parallel::mcmapply(process_covariates, list(X, X_noise, X_scale), list(y,y,y), SIMPLIFY = F)
  names(covariates) = c("X", "X_noise", "X_scale")
  if(any(!covariates$X$no_NA_pattern))message(
    "There are some NAs in X. Whenever possible, 
  e.g. in the case of an intercept or basis functions, 
  put values of X for all times and all places, 
  this will improve greatly the behavior of the associated regression coefficients.
  ")
  ################
  # useful stuff #
  ################
  useful_stuff = list()
  useful_stuff$n_loc = nrow(locs)
  useful_stuff$time_depth = time_depth
  useful_stuff$buffer_depth  = (5*time_depth)*(time_depth!=1)
  useful_stuff$y_split = apply(y, c(1, 3), c, simplify = F) # split  by time and loc for density computation
  useful_stuff$y_na_killed = y; useful_stuff$y_na_killed[is.na(useful_stuff$y_na_killed)] = 0
  useful_stuff$non_na_y = apply(y, c(1, 3), function(x)which(!is.na(x)), simplify = F) 
  useful_stuff$n_var_y = dim(y)[2] # number of variables in y
  useful_stuff$n_var_X = dim(X)[2] # ...
  useful_stuff$n_var_X_noise = dim(X_noise)[2] # ..
  useful_stuff$n_var_X_scale = dim(X_scale)[2] # ..
  useful_stuff$n_time_periods = dim(y)[3] # ..
  useful_stuff$n_field = length(Vecchia_approx_DAG$field_position$location_idx) # number of loc-var pairs in the latent field
  useful_stuff$non_na_count = apply(y, 2, function(x)sum(!is.na(x))) # number of non na obs per space-time position
  # possible configurations of na patterns in y
  useful_stuff$y_NA_possibilities = as.matrix(expand.grid(lapply(seq(dim(y)[2]), function(i)c(T,F))))#[-2^useful_stuff$n_var_y,]
  useful_stuff$y_NA_possibilities_match = as.matrix(array_matrix_mult(is.na(y), matrix(2^seq(dim(y)[2]-1, 0)))[,1,])+1
  useful_stuff$lower_tri_idx = get_lower_tri_idx_DAG(Vecchia_approx_DAG$DAG) # lower triangular idx for Vecchia approx
  useful_stuff$var_idx = get_var_idx(Vecchia_approx_DAG$DAG, Vecchia_approx_DAG$field_position$var_idx) # lower triangular idx for Veccchia approx
  # var-loc couples of y where there is at least one observation in all time periods
  useful_stuff$y_at_least_one_obs = apply(y, c(1, 2), function(x)!all(is.na(x))); useful_stuff$y_at_least_one_obs = split(useful_stuff$y_at_least_one_obs, row(useful_stuff$y_at_least_one_obs)) # loc - var pairs with at least one obs
  useful_stuff$n_field_per_site = lapply(useful_stuff$y_at_least_one_obs, sum)# number of latent field variables simulated per spatial site
  # position of current observation in the var-loc couples of y where there is at least one observation in all time periods
  useful_stuff$position_in_y_at_least_one_obs = 
    apply(y, 3, function(y_slice)
    {
      mapply(
        FUN = function(y_obs, y_at_least_one_obs)return(which(!is.na(y_obs[y_at_least_one_obs]))),
        y_obs = split(y_slice, row(y_slice)),
        y_at_least_one_obs = useful_stuff$y_at_least_one_obs, 
        SIMPLIFY = F
      )
    },
    simplify = F
    )
  
  ######################
  # Hierarchical model #
  ######################
  hierarchical_model = list()
  smoothness_prior = NULL
  if(is.null(smoothness_prior)){
    message("putting a wide prior (.5 <-> 2) on smoothness parameters")
    hierarchical_model$smoothness_prior = cbind("min" = rep(.5, useful_stuff$n_var_y), "max" = rep(2.5, useful_stuff$n_var_y))
  }
  range_prior = NULL
  if(is.null(range_prior)){
    message("putting a wide prior (0 <-> spatial domain radius) on range parameters")
    hierarchical_model$range_prior = cbind("min" = rep(0, useful_stuff$n_var_y), "max" = rep(max(fields::rdist(matrix(apply(locs, 2, mean),1), locs)), useful_stuff$n_var_y))
  }
  scale_beta_prior = NULL
  if(is.null(scale_beta_prior)){
    message("putting a loose Normal prior on scale parameters")
    hierarchical_model$scale_beta_prior = list(
      mean = rep(0, useful_stuff$n_var_X_scale * useful_stuff$n_var_y), 
      sd   = diag(.01, useful_stuff$n_var_X_scale * useful_stuff$n_var_y)
      )
  }
  noise_beta_prior = NULL
  if(is.null(noise_beta_prior)){
    message("putting a loose Normal prior on noise parameters")
    hierarchical_model$noise_beta_prior = list(
      mean = rep(0, useful_stuff$n_var_X_noise    * useful_stuff$n_var_y*(useful_stuff$n_var_y+1)/2), 
      sd   = diag(.01, useful_stuff$n_var_X_noise * useful_stuff$n_var_y*(useful_stuff$n_var_y+1)/2)
      )
  }
  ################
  # chain states # 
  ################
  chains = list()
  for(i in seq(n_chains))
  {
    chains[[i]] = list()
    # model parameters
    chains[[i]]$params = list()
    chains[[i]]$params$noise_beta = matrix(0, useful_stuff$n_var_X_noise, useful_stuff$n_var_y*(useful_stuff$n_var_y+1)/2)#; chains[[i]]$params$noise_beta[] = rnorm(length(chains[[i]]$params$noise_beta))
    chains[[i]]$params$scale_beta = matrix(0, useful_stuff$n_var_X_scale, useful_stuff$n_var_y)                           #; chains[[i]]$params$scale_beta[] = rnorm(length(chains[[i]]$params$scale_beta))
    chains[[i]]$params$beta       = matrix(0, useful_stuff$n_var_X      , useful_stuff$n_var_y)                           ; chains[[i]]$params$beta      [] = rnorm(length(chains[[i]]$params$beta      ))
    chains[[i]]$params$range      = matrix(runif(useful_stuff$n_var_y, min = hierarchical_model$range_prior[,1],      max = hierarchical_model$range_prior[,2]))
    chains[[i]]$params$smoothness = matrix(runif(useful_stuff$n_var_y, min = hierarchical_model$smoothness_prior[,1], max = hierarchical_model$smoothness_prior[,2]))
    chains[[i]]$params$field = matrix(0, useful_stuff$n_time_periods, useful_stuff$n_field)
    rho = GpGp::exponential_isotropic(c(1, 1, 0), matrix(rnorm(2*useful_stuff$n_var_y), useful_stuff$n_var_y))
    chains[[i]]$params$rho_vec = rho[lower.tri(rho, diag = F)];  remove(rho)
    chains[[i]]$params$a2_vec = 1/runif(nrow(hierarchical_model$range_prior), min = hierarchical_model$range_prior[,1], max = hierarchical_model$range_prior[,2])
    chains[[i]]$params$nu_vec = .5 + 2*runif(useful_stuff$n_var_y)
    if(time_depth>1)
    {
      chains[[i]]$params$alpha = .01 * runif(1)
      chains[[i]]$params$a_scal = runif(1) 
      chains[[i]]$params$b_scal = runif(1) 
      chains[[i]]$params$cc = 1*runif(1) 
      chains[[i]]$params$delta = runif(1) 
      chains[[i]]$params$r = 1 * runif(1) 
      chains[[i]]$params$A_vec = runif(useful_stuff$n_var_y)
      chains[[i]]$params$lambda = runif(1) 
    }
    
    chains[[i]]$stuff = list()
    chains[[i]]$stuff$noise_info = get_noise_info(
      X_noise_list = covariates$X_noise, 
      noise_beta = chains[[i]]$params$noise_beta, 
      y_NA_possibilities_match = useful_stuff$y_NA_possibilities_match, 
      y_NA_possibilities = useful_stuff$y_NA_possibilities, 
      y = y)
    
    chains[[i]]$stuff$vecchia = 
      vecchia_block_approx( 
    Vecchia_approx_DAG = Vecchia_approx_DAG, locs = locs, lower_tri_idx_DAG = useful_stuff$lower_tri_idx, 
    var_idx = useful_stuff$var_idx, time_depth = useful_stuff$time_depth, #does not depend on params
    rho_vec = chains[[i]]$params$rho_vec, a =  chains[[i]]$params$a_scal,
    b = chains[[i]]$params$b_scal, cc = chains[[i]]$params$cc, delta = chains[[i]]$params$delta, 
    lambda = chains[[i]]$params$lambda, r = chains[[i]]$params$r, 
    A_vec = chains[[i]]$params$A_vec, nu_vec = chains[[i]]$params$nu_vec, a2_vec = chains[[i]]$params$a2_vec
      )
  }
  names(chains) = paste("chain", seq(n_chains), sep = "_")
  
  list("chains" = chains, "useful_stuff" = useful_stuff, "covariates"= covariates, "Vecchia_approx_DAG" = Vecchia_approx_DAG, "y" = y, "locs" = locs, "hierarchical_model" = hierarchical_model)
}
