
source("MultiNNGP/R/GMA.R")
source("MultiNNGP/R/multivariate_NN.R")
source("MultiNNGP/R/vecchia.R")
source("MultiNNGP/R/initialize.R")
source("MultiNNGP/R/latent_field_sampling.R")



# simulating a latent field

n_loc = 1000
n_var = 3
n_time = 1

rho_vec = rep(.8, n_var*(n_var-1)/2)
nu_vec = c(.5, 1.3, 2.5)
a2_vec = c(10,100,1000)
locs = cbind(1*runif(n_loc), 1)


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


y_true = array(dim = c(n_loc, n_var, 1))
y_true[] = t(matrix(vecchia_blocks_solve(vecchia_blocks = vecchia_blocks, x = rnorm(n_loc*n_var), time_depth = 1),3))



plot(rep(locs[,1],3), y_true, col = rep(seq(3), each = nrow(locs)), pch  = 16, cex = .5)

y = y_true + .3*rnorm(length(y_true))
y[sample(x = seq(length(y)), size = 300)] = NA

plot(rep(locs[,1],3), y, col = rep(seq(3), each = nrow(locs)), pch  = 16, cex = .5)

mcmc_nngp_list = multivariate_NNGP_initialize(
  y = y, locs = locs, X = X, X_noise = X_noise, X_scale = X_scale, Vecchia_approx_DAG = Vecchia_approx_DAG, n_chains = 2)





# getting stuff
vecchia_blocks = mcmc_nngp_list$chains[[1]]$stuff$vecchia
precision_blocks = get_precision_blocks(vecchia_blocks)
noise_info = mcmc_nngp_list$chains[[1]]$stuff$noise_info
X_noise_list = mcmc_nngp_list$covariates$X_noise
diag_precision_block = precision_blocks[[1]]
useful_stuff = mcmc_nngp_list$useful_stuff

latent_field = mcmc_nngp_list$chains$chain_1$params$field

loc_subsetting = (locs[,1]>.5)
plot(locs[,1], loc_subsetting)

# subsetting locations
for(i in unique(loc_subsetting))
{
  loc_subset_idx = which(loc_subsetting==i)
  field_subset_idx = which(mcmc_nngp_list$Vecchia_approx_DAG$field_position$location_idx  %in% loc_subset_idx)
  a_posteriori_diag_precision_chols = get_a_posteriori_diag_precision_chols(
      diag_precision_block, 
      noise_info, 
      X_noise_list, 
      useful_stuff, 
      loc_subset_idx, 
      Vecchia_approx_DAG
  )
  a_posteriori_diag_precision = get_a_posteriori_diag_precision(
      diag_precision_block, 
      noise_info, 
      X_noise_list, 
      useful_stuff, 
      seq(1000), 
      Vecchia_approx_DAG
  )
  field_obs_precision = get_field_obs_precision(
    noise_info, 
    X_noise_list, 
    useful_stuff, 
    loc_subset_idx, 
    Vecchia_approx_DAG, 
    y = mcmc_nngp_list$y
  )
  
  latent_field[field_subset_idx] = 
  posterior_precision_solve_chol(
  a_posteriori_diag_precision_chols = a_posteriori_diag_precision_chols, 
  below_diag_precision_blocks = NULL, 
  vecchia_blocks = vecchia_blocks, 
  loc_subset_idx = loc_subset_idx, 
  x = rnorm(field_subset_idx) + 
    posterior_precision_solve_t_chol(
    a_posteriori_diag_precision_chols = a_posteriori_diag_precision_chols, 
    below_diag_precision_blocks = NULL, 
    vecchia_blocks = vecchia_blocks, 
    loc_subset_idx = loc_subset_idx, 
    x = a_posteriori_diag_precision[[1]][field_subset_idx,-field_subset_idx]%*%latent_field[-field_subset_idx]
    ) + 
    
  
  )
  
  
  
  plot(rep(locs[,1],3), y, col = rep(seq(3), each = nrow(locs)), pch  = 16, cex = .5)
  points(mcmc_nngp_list$locs[Vecchia_approx_DAG$field_position$location_idx, 1], latent_field,  col = Vecchia_approx_DAG$field_position$var_idx)
  
}
















