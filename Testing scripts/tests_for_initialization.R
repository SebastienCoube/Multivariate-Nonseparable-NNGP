source("MultiNNGP/R/GMA.R")
source("MultiNNGP/R/multivariate_NN.R")
source("MultiNNGP/R/vecchia.R")
source("MultiNNGP/R/initialize.R")



n_loc = 1000
n_var = 5
n_time = 100

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
Vecchia_approx_DAG = make_simple_Vecchia_approx_DAG(
  y = y, locs = locs
)
mcmc_nngp_list = multivariate_NNGP_initialize(
  y = y, locs = locs, X = X, X_noise = X_noise, X_scale = X_scale, Vecchia_approx_DAG = Vecchia_approx_DAG, n_chains = 2)

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
Vecchia_approx_DAG = make_simple_Vecchia_approx_DAG(
  y = y, locs = locs
)
mcmc_nngp_list = multivariate_NNGP_initialize(
  y = y, locs = locs, X = X, X_noise = X_noise, X_scale = X_scale, Vecchia_approx_DAG = Vecchia_approx_DAG, n_chains = 2)

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
Vecchia_approx_DAG = make_simple_Vecchia_approx_DAG(
  y = y, locs = locs
)
mcmc_nngp_list = multivariate_NNGP_initialize(
  y = y, locs = locs, X = X, X_noise = X_noise, X_scale = X_scale, Vecchia_approx_DAG = Vecchia_approx_DAG, n_chains = 2)

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
Vecchia_approx_DAG = make_simple_Vecchia_approx_DAG(
  y = y, locs = locs
)
mcmc_nngp_list = multivariate_NNGP_initialize(
  y = y, locs = locs, X = X, X_noise = X_noise, X_scale = X_scale, Vecchia_approx_DAG = Vecchia_approx_DAG, n_chains = 2)

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
Vecchia_approx_DAG = make_simple_Vecchia_approx_DAG(
  y = y, locs = locs
)
mcmc_nngp_list = multivariate_NNGP_initialize(
  y = y, locs = locs, X = X, X_noise = X_noise, X_scale = X_scale, Vecchia_approx_DAG = Vecchia_approx_DAG, n_chains = 2)



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
Vecchia_approx_DAG = make_simple_Vecchia_approx_DAG(
  y = y, locs = locs
)
mcmc_nngp_list = multivariate_NNGP_initialize(
  y = y, locs = locs, X = X, X_noise = X_noise, X_scale = X_scale, Vecchia_approx_DAG = Vecchia_approx_DAG, n_chains = 2)



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
Vecchia_approx_DAG = make_simple_Vecchia_approx_DAG(
  y = y, locs = locs
)
mcmc_nngp_list = multivariate_NNGP_initialize(
  y = y, locs = locs, X = X, X_noise = X_noise, X_scale = X_scale, Vecchia_approx_DAG = Vecchia_approx_DAG, n_chains = 2)

#################################################################################
# test : some loc-var couples never observed in Y, only an intercept in X_noise #
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
X_noise[,,seq(50)] = NA
X_noise = X_noise[,1,,drop = F]
X_scale = X
locs = matrix(runif(2*n_loc), n_loc)
Vecchia_approx_DAG = make_simple_Vecchia_approx_DAG(
  y = y, locs = locs
)
mcmc_nngp_list = multivariate_NNGP_initialize(
  y = y, locs = locs, X = X, X_noise = X_noise, X_scale = X_scale, Vecchia_approx_DAG = Vecchia_approx_DAG, n_chains = 2)

#################################################################################
# test : some loc-var couples never observed in Y, only an intercept in X_noise #
#################################################################################
y = array(dim = c(n_loc, n_var, n_time))
dimnames(y) = list(
  paste("loc", seq(n_loc), sep = "_"), 
  paste("var", seq(n_var), sep = "_"), 
  paste("time", seq(n_time), sep = "_")
)
y[] = rnorm(length(y))
y[,,seq(50)] = NA

y[,,seq(61, 70)] = NA

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
X_noise = X_noise[,1,,drop = F]
X_scale = X
locs = matrix(runif(2*n_loc), n_loc)
Vecchia_approx_DAG = make_simple_Vecchia_approx_DAG(
  y = y, locs = locs
)
mcmc_nngp_list = multivariate_NNGP_initialize(
  y = y, locs = locs, X = X, X_noise = X_noise, X_scale = X_scale, Vecchia_approx_DAG = Vecchia_approx_DAG, n_chains = 2)



#####################
# test : only space #
#####################
y = array(dim = c(n_loc, n_var, 1))
dimnames(y) = list(
  paste("loc", seq(n_loc), sep = "_"), 
  paste("var", seq(n_var), sep = "_"), 
  paste("time", seq(1), sep = "_")
)
y[] = rnorm(length(y))
y[cbind(
  1+floor(dim(y)[1]*runif(1000)),
  1+floor(dim(y)[2]*runif(1000)),
  1+floor(dim(y)[3]*runif(1000))
)] = NA
X = array(dim = c(n_loc, 10, 1))
dimnames(X) = list(
  paste("loc", seq(n_loc), sep = "_"), 
  paste("covariate", seq(dim(X)[2]), sep = "_"), 
  paste("time", seq(1), sep = "_")
)
X[] = rnorm(length(X))
X[,1,] = 1
X_noise = X
X_scale = X
locs = matrix(runif(2*n_loc), n_loc)
locs = cbind(1000*runif(1*n_loc), 0*runif(1*n_loc))
Vecchia_approx_DAG = make_simple_Vecchia_approx_DAG(
  y = y, locs = locs
)
mcmc_nngp_list = multivariate_NNGP_initialize(
  y = y, locs = locs, X = X, X_noise = X_noise, X_scale = X_scale, Vecchia_approx_DAG = Vecchia_approx_DAG, n_chains = 2)


i=150
plot(locs[Vecchia_approx_DAG$field_position$location_idx[seq(i)],])
points(locs[Vecchia_approx_DAG$field_position$location_idx[Vecchia_approx_DAG$DAG$children[[(i)]]],,drop=F], pch=3)
points(locs[Vecchia_approx_DAG$field_position$location_idx[Vecchia_approx_DAG$DAG$parents_same_time[[(i)]]],,drop=F], pch=4)
