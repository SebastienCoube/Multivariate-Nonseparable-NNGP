### source("MultiNNGP/R/GMA.R")
### source("MultiNNGP/R/multivariate_NN.R")
### source("MultiNNGP/R/vecchia.R")
 my_image = function(m)image(t(m)[,nrow(m):1])


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
  # adding buffer at the beginning
  if(dim(y)[3]>1){
    X =       abind::abind(array(123456 , c(dim(X)      [c(1,2)], 5*time_depth)), X,       along = 3)
    X_scale = abind::abind(array(123456 , c(dim(X_scale)[c(1,2)], 5*time_depth)), X_scale, along = 3)
    X_noise = abind::abind(array(NA,      c(dim(X_noise)[c(1,2)], 5*time_depth)), X_noise, along = 3)
    y =       abind::abind(array(NA ,     c(dim(y)      [c(1,2)], 5*time_depth)), y,       along = 3)
  }
  # adding buffer at the end
  if(dim(y)[3]>1){
    X =       abind::abind(X,       array(123456 , c(dim(X)      [c(1,2)], 5*time_depth)), along = 3)
    X_scale = abind::abind(X_scale, array(123456 , c(dim(X_scale)[c(1,2)], 5*time_depth)), along = 3)
    X_noise = abind::abind(X_noise, array(NA,      c(dim(X_noise)[c(1,2)], 5*time_depth)), along = 3)
    y =       abind::abind(y,       array(NA ,     c(dim(y)      [c(1,2)], 5*time_depth)), along = 3)
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
###  # var-loc couples of y where there is at least one observation in all time periods
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
  useful_stuff$y_loc_var_format = t(mapply(function(i, j)y[i,j,], i = Vecchia_approx_DAG$field_position$location_idx, j = Vecchia_approx_DAG$field_position$var_idx))
  useful_stuff$y_loc_var_format_no_NA = useful_stuff$y_loc_var_format; useful_stuff$y_loc_var_format_no_NA[is.na(useful_stuff$y_loc_var_format_no_NA)]=0
  useful_stuff$noise_precision_i = lapply(
    seq(useful_stuff$n_time_periods), function(i_time)
      unlist(lapply(
        seq(useful_stuff$n_loc), function(i_loc){
          tatato = Vecchia_approx_DAG$field_position$loc_match[[i_loc]][useful_stuff$position_in_y_at_least_one_obs[[i_time]][[i_loc]]]
          rep(tatato, length(tatato))
        }
      ))
  )
  useful_stuff$noise_precision_j = lapply(
    seq(useful_stuff$n_time_periods), function(i_time)
      unlist(lapply(
        seq(useful_stuff$n_loc), function(i_loc){
          tatato = Vecchia_approx_DAG$field_position$loc_match[[i_loc]][useful_stuff$position_in_y_at_least_one_obs[[i_time]][[i_loc]]]
          rep(tatato, each = length(tatato))
        }
      ))
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
  
  list("useful_stuff" = useful_stuff, "covariates"= covariates, "Vecchia_approx_DAG" = Vecchia_approx_DAG, "y" = y, "locs" = locs, "hierarchical_model" = hierarchical_model)
}



make_one_chain = function(mcmc_nngp_list)
{
  ################
  # chain states # 
  ################
  chain = list()
  # model parameters
  chain$params = list()
  chain$params$noise_beta = matrix(0, useful_stuff$n_var_X_noise, useful_stuff$n_var_y*(useful_stuff$n_var_y+1)/2)#; chain$params$noise_beta[] = rnorm(length(chain$params$noise_beta))
  chain$params$scale_beta = matrix(0, useful_stuff$n_var_X_scale, useful_stuff$n_var_y)                           #; chain$params$scale_beta[] = rnorm(length(chain$params$scale_beta))
  chain$params$beta       = matrix(0, useful_stuff$n_var_X      , useful_stuff$n_var_y)                           ; chain$params$beta      [] = rnorm(length(chain$params$beta      ))
  chain$params$range      = matrix(runif(useful_stuff$n_var_y, min = hierarchical_model$range_prior[,1],      max = hierarchical_model$range_prior[,2]))
  chain$params$smoothness = matrix(runif(useful_stuff$n_var_y, min = hierarchical_model$smoothness_prior[,1], max = hierarchical_model$smoothness_prior[,2]))
  chain$params$field = matrix(0, useful_stuff$n_field, useful_stuff$n_time_periods)
  rho = GpGp::exponential_isotropic(c(1, 1, 0), matrix(rnorm(2*useful_stuff$n_var_y), useful_stuff$n_var_y))
  chain$params$rho_vec = rho[lower.tri(rho, diag = F)];  remove(rho)
  chain$params$a2_vec = 1/runif(nrow(hierarchical_model$range_prior), min = hierarchical_model$range_prior[,1], max = hierarchical_model$range_prior[,2])
  chain$params$nu_vec = .5 + 2*runif(useful_stuff$n_var_y)
  if(time_depth>1)
  {
    chain$params$alpha = .01 * runif(1)
    chain$params$a_scal = runif(1) 
    chain$params$b_scal = runif(1) 
    chain$params$cc = 1*runif(1) 
    chain$params$delta = runif(1) 
    chain$params$r = 1 * runif(1) 
    chain$params$A_vec = runif(useful_stuff$n_var_y)
    chain$params$lambda = runif(1) 
  }
  
  chain$stuff = list()
  chain$stuff$noise_info = get_noise_info(
    X_noise_list = covariates$X_noise, 
    noise_beta = chain$params$noise_beta, 
    y_NA_possibilities_match = useful_stuff$y_NA_possibilities_match, 
    y_NA_possibilities = useful_stuff$y_NA_possibilities, 
    y = y)
  
  chain$stuff$vecchia = 
    vecchia_block_approx( 
      Vecchia_approx_DAG = Vecchia_approx_DAG, locs = locs, lower_tri_idx_DAG = useful_stuff$lower_tri_idx, 
      var_idx = useful_stuff$var_idx, time_depth = useful_stuff$time_depth, #does not depend on params
      rho_vec = chain$params$rho_vec, a =  chain$params$a_scal,
      b = chain$params$b_scal, cc = chain$params$cc, delta = chain$params$delta, 
      lambda = chain$params$lambda, r = chain$params$r, 
      A_vec = chain$params$A_vec, nu_vec = chain$params$nu_vec, a2_vec = chain$params$a2_vec
    )
  return(chain)
}