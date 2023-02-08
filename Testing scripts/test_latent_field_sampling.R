
source("MultiNNGP/R/GMA.R")
source("MultiNNGP/R/multivariate_NN.R")
source("MultiNNGP/R/vecchia.R")
source("MultiNNGP/R/initialize.R")
source("MultiNNGP/R/latent_field_sampling.R")



# simulating a latent field

n_loc = 1000
n_var = 3
n_time = 1

rho_vec = rep(.95, n_var*(n_var-1)/2)
nu_vec = c(.5, 1.3, 2.5)
a2_vec = c(10,20,30)
locs = cbind(10*runif(n_loc), 1)


X = array(dim = c(n_loc, 1, 1))
dimnames(X) = list(
  paste("loc", seq(n_loc), sep = "_"), 
  paste("covariate", seq(dim(X)[2]), sep = "_"), 
  paste("time", seq(1), sep = "_")
)
X[] = rnorm(length(X))
X[,1,] = 1
X_noise = X
X_scale = X
Vecchia_approx_DAG = make_simple_Vecchia_approx_DAG(
  y = array(1, dim = c(n_loc, n_var, 1)), locs = locs
)

lower_tri_idx_DAG = get_lower_tri_idx_DAG(Vecchia_approx_DAG$DAG)
var_idx = get_var_idx(Vecchia_approx_DAG$DAG, Vecchia_approx_DAG$field_position$var_idx)
vecchia_blocks = vecchia_block_approx(
    Vecchia_approx_DAG = Vecchia_approx_DAG, locs = locs, 
    lower_tri_idx_DAG = lower_tri_idx_DAG, 
    var_idx= var_idx, time_depth = 1, #does not depend on params
    rho_vec = rho_vec, a = NULL, b = NULL, cc = NULL, delta = NULL, lambda = NULL, r = NULL, 
    A_vec = NULL, nu_vec = nu_vec, a2_vec = a2_vec
)


y = array(dim = c(n_loc, n_var, 1))
y[] = t(matrix(vecchia_blocks_solve(vecchia_blocks = vecchia_blocks, x = rnorm(n_loc*n_var), time_depth = 1),3))

plot(locs[,1], y[,1,])
points(locs[,1], y[,2,], col=2)
points(locs[,1], y[,3,], col=3)



mcmc_nngp_list = multivariate_NNGP_initialize(
  y = y, locs = locs, X = X, X_noise = X_noise, X_scale = X_scale, Vecchia_approx_DAG = Vecchia_approx_DAG, n_chains = 2)





# getting stuff
vecchia_blocks = mcmc_nngp_list$chains[[1]]$stuff$vecchia
precision_blocks = get_precision_blocks(vecchia_blocks)
noise_info = mcmc_nngp_list$chains[[1]]$stuff$noise_info
X_noise_list = mcmc_nngp_list$covariates$X_noise
diag_precision_block = precision_blocks[[1]]
useful_stuff = mcmc_nngp_list$useful_stuff

# subsetting locations
loc_subset_idx = seq(useful_stuff$n_loc/2)








