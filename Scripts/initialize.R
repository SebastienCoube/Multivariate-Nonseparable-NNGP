problem with noise info length in mcmc nngp list
# generating data 
n_loc = 500
n_time = 200
n_var = 5

# full data
y = array(dim = c(n_loc, n_var, n_time))
dimnames(y) = list(
  paste("loc", seq(n_loc), sep = "_"), 
  paste("var", seq(n_var), sep = "_"), 
  paste("time", seq(n_time), sep = "_")
  )
y[] = rnorm(length(y))

X = array(dim = c(n_loc, 10, n_time))
dimnames(X) = list(
  paste("loc", seq(n_loc), sep = "_"), 
  paste("covariate", seq(dim(X)[2]), sep = "_"), 
  paste("time", seq(n_time), sep = "_")
  )
X[] = rnorm(length(X))
X[,1,] = 1
X_noise = X
X_scale = X

locs = matrix(runif(2*n_loc), n_loc)

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
  res = list()
  expmat_and_derivatives = get_expmat_and_derivatives(mat_coordinates = mat_coordinates)
  res$tau_precision = get_tau_precision(expmat_and_derivatives$expmat, y)
  res$tau_precision_deriv = lapply(expmat_and_derivatives$expmat_derivative, function(x) res$tau_precision %*% x[!is.na(y), !is.na(y)] %*% res$tau_precision)
  return(res)
}


# get all noise precision matrices and their derivatives following NA pattern of y 
get_noise_info = function(X_noise_list, noise_beta, y, y_NA_possibilities)
{
  res = list()
  mat_coordinates_array = array_matrix_mult(X_noise_list$X_killed_NA, noise_beta)
  # nonstationary case
  if(!X_noise_list$is_only_intercept) expmat_and_derivatives = parallel::mcmapply(
    FUN = get_tau_info, 
    mat_coordinates = lapply(X_noise_list$no_NA_idx_split, function(x)mat_coordinates_array[x[1],,x[2]]),
    y = lapply(X_noise_list$no_NA_idx_split, function(x)y[x[1],,x[2]]), 
    SIMPLIFY = F, mc.cores = parallel::detectCores()-1
  )
  # stationary case
  if(X_noise_list$is_only_intercept)
  {
    tau_info_cases = lapply(
      split(y_NA_possibilities, row(y_NA_possibilities)), 
      function(x)get_tau_info(mat_coordinates = rep(1, nrow(noise_beta))%*%noise_beta, y = x)
    )
    expmat_and_derivatives = lapply(
      X_noise_list$y_NA_pattern, 
      function(x)tau_info_cases[[x]])
  }
    
  return(expmat_and_derivatives)
  #res$tau_precision = get_tau_precision(expmat_and_derivatives$expmat, y)
  #res$tau_precision_deriv = lapply(expmat_and_derivatives$expmat_derivative, function(x) res$tau_precision %*% x[!is.na(y), !is.na(y)] %*% res$tau_precision)
  #res
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
  res$no_NA_idx = cbind(row(res$no_NA_pattern)[res$no_NA_pattern], col(res$no_NA_pattern)[res$no_NA_pattern]) # dim 1 and 3 of complete observations 
  res$no_NA_idx_split = split(res$no_NA_idx, row(res$no_NA_idx))  # dim 1 and 3 of complete observations, split into list
  res$is_only_intercept = (dim(X)[2]==1)&all(na.omit(X)==1) # bool to check if X is only an intercept
  
  # used in intercept only. tells which y are observed at thee indices of res$no_NA_idx_split
  y_NA_possibilities = as.matrix(expand.grid(lapply(seq(dim(y)[2]), function(i)c(1,NA))))
  y_NA_possibilities = y_NA_possibilities[-nrow(y_NA_possibilities),]
  y_NA_possibilities = split(y_NA_possibilities, row(y_NA_possibilities))
  res$y_NA_pattern = match(
    lapply(res$no_NA_idx_split, function(x)is.na(unname(y[x[1],,x[2]]))[seq(dim(y)[2])]),
    lapply(y_NA_possibilities, is.na)
  )
  res
}


make_simple_DAG = 
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
  NNarray_pevious_times = find_unordered_nn_multi(
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
        "var_idx" = var_tag[at_least_one_obs]
      )
    )
  )
  }

DAG = make_simple_DAG(
  y = y, locs = locs
)



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
    DAG, 
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
  if(any(dim(y)[c(1,3)]!=dim(X)[c(1,3)]))stop("dimensions 1 (spatial location), 3 (time) and 4(observations) of y and X should be the same")
  # checking X
  if(any((apply(y, c(1, 3), anyNA)==F)&(apply(X, c(1, 3), anyNA)==T)))stop("X should have no NAs where y has no NAs (but it is possible, and even recommended, for X to have observations even where y has NAs)")
  if(any((apply(y, c(1, 3), anyNA)==F)&(apply(X_noise, c(1, 3), anyNA)==T)))stop("X_noise should have no NAs where y has no NAs")
  if(any(is.na(X_scale)))stop("X_scale can have no NA")
  #########################
  # Processing covariates #
  #########################
  covariates = list(X = process_covariates(X, y), X_noise = process_covariates(X_noise, y), X_scale = process_covariates(X_scale, y))
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
  useful_stuff$y_split = apply(y, c(1, 3), c, simplify = F) # split  by time and loc for density computation
  useful_stuff$non_na_y = apply(y, c(1, 3), function(x)which(!is.na(x)), simplify = F) 
  useful_stuff$n_var_y = dim(y)[2] # number of variables in y
  useful_stuff$n_var_X = dim(X)[2] # ...
  useful_stuff$n_var_X_noise = dim(X_noise)[2] # ..
  useful_stuff$n_var_X_scale = dim(X_scale)[2] # ..
  useful_stuff$n_time_periods = dim(y)[3] # ..
  useful_stuff$n_field = length(DAG$field_position$location_idx) # number of loc-var pairs in the latent field
  useful_stuff$non_na_count = apply(y, 2, function(x)sum(!is.na(x))) # number of non na obs per space-time position
  # possible configurations of na patterns in y
  useful_stuff$y_NA_possibilities = as.matrix(expand.grid(lapply(seq(dim(y)[2]), function(i)c(1,NA))))[-2^useful_stuff$n_var_y,]
  useful_stuff$lower_tri_idx = get_lower_tri_idx_DAG(DAG$DAG) # lower triangular idx for Veccchia approx
  useful_stuff$var_idx = get_var_idx(DAG$DAG, DAG$field_position$var_idx) # lower triangular idx for Veccchia approx
  
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
    
    chains[[i]]$stuff = list()
    chains[[i]]$stuff$noise_info = get_noise_info(
      X_noise_list = covariates$X_noise, 
      noise_beta = chains[[i]]$params$noise_beta, 
      y = y, y_NA_possibilities = useful_stuff$y_na_possibilities)
    chains[[i]]$stuff$vecchia = 
  }
  
}


########################
# test : some NAs in Y #
########################
y = array(dim = c(n_loc, n_var, n_time))
dimnames(y) = list(
  paste("loc", seq(n_loc), sep = "_"), 
  paste("var", seq(n_var), sep = "_"), 
  paste("time", seq(n_time), sep = "_")
)
y[] = rnorm(length(y))
y[cbind(
  1+floor(dim(y)[1]*runif(10000)),
  1+floor(dim(y)[2]*runif(10000)),
  1+floor(dim(y)[3]*runif(10000))
)] = NA
X = array(dim = c(n_loc, 10, n_time))
dimnames(X) = list(
  paste("loc", seq(n_loc), sep = "_"), 
  paste("covariate", seq(dim(X)[2]), sep = "_"), 
  paste("time", seq(n_time), sep = "_")
)
X[] = rnorm(length(X))
X[,1,] = 1
X_noise = X
X_scale = X
locs = matrix(runif(2*n_loc), n_loc)
DAG = make_simple_DAG(
  y = y, locs = locs
)
mcmc_nngp_list = multivariate_NNGP_initialize(
  y = y, locs = locs, X = X, X_noise = X_noise, X_scale = X_scale, DAG = DAG, n_chains = 2)

###################################################
# test : some loc-var couples never observed in Y #
###################################################
y = array(dim = c(n_loc, n_var, n_time))
dimnames(y) = list(
  paste("loc", seq(n_loc), sep = "_"), 
  paste("var", seq(n_var), sep = "_"), 
  paste("time", seq(n_time), sep = "_")
)
y[] = rnorm(length(y))
y[cbind(
  1+floor(dim(y)[1]*runif(10000)),
  1+floor(dim(y)[2]*runif(10000)),
  1+floor(dim(y)[3]*runif(10000))
)] = NA
y[seq(10)    ,1,] = NA
y[seq(11, 20),2,] = NA
y[seq(21, 30),3,] = NA
y[seq(31, 40),4,] = NA
y[seq(41, 50),5,] = NA

X = array(dim = c(n_loc, 10, n_time))
dimnames(X) = list(
  paste("loc", seq(n_loc), sep = "_"), 
  paste("covariate", seq(dim(X)[2]), sep = "_"), 
  paste("time", seq(n_time), sep = "_")
)
X[] = rnorm(length(X))
X[,1,] = 1
X_noise = X
X_scale = X
locs = matrix(runif(2*n_loc), n_loc)
DAG = make_simple_DAG(
  y = y, locs = locs
)
mcmc_nngp_list = multivariate_NNGP_initialize(
  y = y, locs = locs, X = X, X_noise = X_noise, X_scale = X_scale, DAG = DAG, n_chains = 2)

############################################################################
# test : some loc-var couples never observed in Y and some hallal NAs in X #
############################################################################
y = array(dim = c(n_loc, n_var, n_time))
dimnames(y) = list(
  paste("loc", seq(n_loc), sep = "_"), 
  paste("var", seq(n_var), sep = "_"), 
  paste("time", seq(n_time), sep = "_")
)
y[] = rnorm(length(y))
y[,,seq(50)] = NA

X = array(dim = c(n_loc, 10, n_time))
dimnames(X) = list(
  paste("loc", seq(n_loc), sep = "_"), 
  paste("covariate", seq(dim(X)[2]), sep = "_"), 
  paste("time", seq(n_time), sep = "_")
)
X[] = rnorm(length(X))
X[,1,] = 1
X_noise = X
X_scale = X
X[,,seq(50)] = NA
locs = matrix(runif(2*n_loc), n_loc)
DAG = make_simple_DAG(
  y = y, locs = locs
)
mcmc_nngp_list = multivariate_NNGP_initialize(
  y = y, locs = locs, X = X, X_noise = X_noise, X_scale = X_scale, DAG = DAG, n_chains = 2)

###########################################################################
# test : some loc-var couples never observed in Y and some haram NAs in X #
###########################################################################
y = array(dim = c(n_loc, n_var, n_time))
dimnames(y) = list(
  paste("loc", seq(n_loc), sep = "_"), 
  paste("var", seq(n_var), sep = "_"), 
  paste("time", seq(n_time), sep = "_")
)
y[] = rnorm(length(y))
y[,,seq(50)] = NA

X = array(dim = c(n_loc, 10, n_time))
dimnames(X) = list(
  paste("loc", seq(n_loc), sep = "_"), 
  paste("covariate", seq(dim(X)[2]), sep = "_"), 
  paste("time", seq(n_time), sep = "_")
)
X[] = rnorm(length(X))
X[,1,] = 1
X_noise = X
X_scale = X
X[,,seq(100)] = NA
locs = matrix(runif(2*n_loc), n_loc)
DAG = make_simple_DAG(
  y = y, locs = locs
)
mcmc_nngp_list = multivariate_NNGP_initialize(
  y = y, locs = locs, X = X, X_noise = X_noise, X_scale = X_scale, DAG = DAG, n_chains = 2)

############################
# test : one NA in X_scale #
############################
y = array(dim = c(n_loc, n_var, n_time))
dimnames(y) = list(
  paste("loc", seq(n_loc), sep = "_"), 
  paste("var", seq(n_var), sep = "_"), 
  paste("time", seq(n_time), sep = "_")
)
y[] = rnorm(length(y))

X = array(dim = c(n_loc, 10, n_time))
dimnames(X) = list(
  paste("loc", seq(n_loc), sep = "_"), 
  paste("covariate", seq(dim(X)[2]), sep = "_"), 
  paste("time", seq(n_time), sep = "_")
)
X[] = rnorm(length(X))
X[,1,] = 1
X_noise = X
X_scale = X
X_scale[1,1,1]=NA
locs = matrix(runif(2*n_loc), n_loc)
DAG = make_simple_DAG(
  y = y, locs = locs
)
mcmc_nngp_list = multivariate_NNGP_initialize(
  y = y, locs = locs, X = X, X_noise = X_noise, X_scale = X_scale, DAG = DAG, n_chains = 2)



##################################################################################
# test : some loc-var couples never observed in Y and some hallal NAs in X_noise #
##################################################################################
y = array(dim = c(n_loc, n_var, n_time))
dimnames(y) = list(
  paste("loc", seq(n_loc), sep = "_"), 
  paste("var", seq(n_var), sep = "_"), 
  paste("time", seq(n_time), sep = "_")
)
y[] = rnorm(length(y))
y[,,seq(50)] = NA

X = array(dim = c(n_loc, 10, n_time))
dimnames(X) = list(
  paste("loc", seq(n_loc), sep = "_"), 
  paste("covariate", seq(dim(X)[2]), sep = "_"), 
  paste("time", seq(n_time), sep = "_")
)
X[] = rnorm(length(X))
X[,1,] = 1
X_noise = X
X_noise[,,seq(50)] = NA
X_scale = X
locs = matrix(runif(2*n_loc), n_loc)
DAG = make_simple_DAG(
  y = y, locs = locs
)
mcmc_nngp_list = multivariate_NNGP_initialize(
  y = y, locs = locs, X = X, X_noise = X_noise, X_scale = X_scale, DAG = DAG, n_chains = 2)



#################################################################################
# test : some loc-var couples never observed in Y and some haram NAs in X_noise #
#################################################################################
y = array(dim = c(n_loc, n_var, n_time))
dimnames(y) = list(
  paste("loc", seq(n_loc), sep = "_"), 
  paste("var", seq(n_var), sep = "_"), 
  paste("time", seq(n_time), sep = "_")
)
y[] = rnorm(length(y))
y[,,seq(50)] = NA

X = array(dim = c(n_loc, 10, n_time))
dimnames(X) = list(
  paste("loc", seq(n_loc), sep = "_"), 
  paste("covariate", seq(dim(X)[2]), sep = "_"), 
  paste("time", seq(n_time), sep = "_")
)
X[] = rnorm(length(X))
X[,1,] = 1
X_noise = X
X_noise[,,seq(100)] = NA
X_scale = X
locs = matrix(runif(2*n_loc), n_loc)
DAG = make_simple_DAG(
  y = y, locs = locs
)
mcmc_nngp_list = multivariate_NNGP_initialize(
  y = y, locs = locs, X = X, X_noise = X_noise, X_scale = X_scale, DAG = DAG, n_chains = 2)


