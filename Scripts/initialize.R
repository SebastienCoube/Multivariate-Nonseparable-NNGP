
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


process_covariates = function(X)
{
  res = list()
  res$X = X 
  res$X_killed_NA = X; res$X_killed_NA[is.na(res$X_killed_NA)] = 0
  res$chol_X = chol(crossprod(apply(res$X_killed_NA, 2, c)))
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
  # checking format
  if(!is.array(y))stop("y should be an array")
  if(!is.array(X))stop("X should be an array")
  if(!is.matrix(locs))stop("locs should be an array")
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
  covariates = list(X = process_covariates(X), X_noise = process_covariates(X_noise), X_scale = process_covariates(X_scale))
  any_missing = apply(X,2, anyNA) 
  if(any(any_missing))message(
    "There are some NAs in X. Whenever possible, 
  e.g. in the case of an intercept or basis functions, 
  put values of X for all times and all places, 
  this will improve greatly the behavior of the associated regression coefficients.
  ")
  covariates$X$complete_obs = seq(dim(X)[2])[!any_missing]
  ################
  # useful stuff #
  ################
  useful_stuff = list()
  useful_stuff$y_split = apply(y, c(1, 3), na.omit, simplify = F) # split  by time and loc for density computation
  useful_stuff$non_na_y = apply(y, c(1, 3), function(x)which(!is.na(x)), simplify = F) 
  useful_stuff$n_var_y = dim(y)[2] # number of variables in y
  useful_stuff$n_var_X = dim(X)[2] # ...
  useful_stuff$n_var_X_noise = dim(X_noise)[2] # ..
  useful_stuff$n_var_X_scale = dim(X_scale)[2] # ..
  useful_stuff$n_time_periods = dim(y)[3] # ..
  useful_stuff$n_field = length(DAG$field_position$location_idx) # number of loc-var pairs in the latent field
  #useful_stuff$field_at_observations_idx = fedwcd
  #tatato = y
  #for(i in seq(length(DAG$field_position$location_idx)))tatato[DAG$field_position$location_idx[i], DAG$field_position$var_idx[i],]=i
  #tatato = apply(tatato, c(1, 3), c, simplify = F)
  #na_array = apply(y, c(1, 3), function(x)which(!is.na(x)), simplify = F)
  #tatata  = tatato
  #for(i in seq(nrow(tatato)))for(j in seq(ncol(tatato)))tatata[i,j][[1]] = tatato[i, j][[1]][na_array[i, j][[1]]]
  #
  #tatata  
  
  ################
  # chain states # 
  ################
  chains = list()
  for(i in seq(n_chains))
  {
    chains[[i]] = list()
    chains[[i]]$params = list()
    chains[[i]]$params$noise_beta = matrix(0, useful_stuff$n_var_X_noise, useful_stuff$n_var_y*(useful_stuff$n_var_y+1)/2); chains[[i]]$params$noise_beta[] = rnorm(length(chains[[i]]$params$noise_beta))
    chains[[i]]$params$scale_beta = matrix(0, useful_stuff$n_var_X_scale, useful_stuff$n_var_y)                           ; chains[[i]]$params$scale_beta[] = rnorm(length(chains[[i]]$params$scale_beta))
    chains[[i]]$params$beta       = matrix(0, useful_stuff$n_var_X      , useful_stuff$n_var_y)                           ; chains[[i]]$params$beta      [] = rnorm(length(chains[[i]]$params$beta      ))
    chains[[i]]$params$field = matrix(0, useful_stuff$n_time_periods, useful_stuff$n_field)
  }
  
  apply(chains[[i]]$params$field, c(1, 2), c, simplify = F)
  
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


