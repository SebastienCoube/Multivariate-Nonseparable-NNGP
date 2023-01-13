
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


process_covariates = function(X)
{
  res = list()
  res$X = X 
  res$X_killed_NA = X; res$X_killed_NA[is.na(res$X_killed_NA)] = 0
  res$chol_X = chol(crossprod(apply(res$X_killed_NA, 2, c)))
  res
}


make_simple_DAG = function(y, locs, 
                           m_whatever_closest_same_time = 10, 
                           m_same_var_same_time = 5, 
                           m_other_vars_same_time = 0, 
                           
                           )
{
  # splitting y and loc in a loc-var format
  var_tag = outer(rep(1, dim(y)[1]), seq(dim(y)[2]))
  loc_idx = outer(seq(dim(y)[1]), rep(1, dim(y)[2]))
  at_least_one_obs = !apply(y, c(1, 2), anyNA)
  locs_ = locs[loc_idx[at_least_one_obs],]
  var_tag = var_tag[at_least_one_obs]
  # reordering
  ordering = GpGp::order_maxmin(locs_)
  #
  NNarray_same_time = find_ordered_nn_multi(locs = locs, 
                                            var_tag = var_tag, 
                                            m_whatever_closest = m_whatever_closest_same_time, 
                                            m_same_var = m_same_var_same_time, 
                                            m_other_vars = m_other_vars_same_time, 
                                            lonlat = F)
  
}
NNarray_pevious_times = find_unordered_nn_multi(locs = locs, 
                                                var_tag = var_tag, 
                                                m_whatever_closest = 10, m_same_var = 5, m_other_vars = 0, 
                                                lonlat = F)
DAG = list(
  "children" = lapply(seq(n_loc), function(x)x), 
  "parents_same_time" = NNarray_same_time, 
  "parents_previous_times" = NNarray_pevious_times 
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
    X_scale
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
  ############################
  # deducting stuff from DAG #
  ############################
  
  
  
}


