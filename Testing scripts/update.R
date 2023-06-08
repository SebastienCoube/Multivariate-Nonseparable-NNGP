
source("MultiNNGP/R/GMA.R")
source("MultiNNGP/R/array_useful_functions.R")
source("MultiNNGP/R/multivariate_NN.R")
source("MultiNNGP/R/vecchia.R")
source("MultiNNGP/R/initialize.R")
source("MultiNNGP/R/latent_field_sampling.R")
source("MultiNNGP/R/basis_functions.R")
source("MultiNNGP/R/noise_variance.R")
#source("MultiNNGP/R/update.R")


# synthetic data set #########

# space time layout
n_loc = 500
n_var = 2
n_time = 50
locs_no_na = cbind(2*runif(n_loc), 1)

# meters
rho_vec = rep(.8, n_var*(n_var-1)/2)
nu_vec = seq(n_var)/n_var+.5
log_range_vec = log(rep(.002, n_var))
A_vec = rep(.5, n_var)
a =       .5
b =       1
delta=    .5
lambda =  .5

cc = 10# .05
r =  10#  .05

# sampling latent field
Vecchia_approx_DAG_no_NA = make_simple_Vecchia_approx_DAG(
  y = array(1, dim = c(n_loc, n_var, n_time)), locs = locs_no_na)
t1 = Sys.time()
vecchia_blocks_no_NA = vecchia_block_approx(
  Vecchia_approx_DAG = Vecchia_approx_DAG_no_NA, 
  locs = locs_no_na, 
  time_depth = 5, #does not depend on params
  rho_vec = rho_vec, a = a, b = b, cc = cc, delta = delta, lambda = lambda, r = r, 
  A_vec = A_vec, nu_vec = nu_vec, log_range_vec = log_range_vec
)
Sys.time()-t1

tatato = vecchia_blocks_solve(vecchia_blocks = vecchia_blocks_no_NA, 
                              x = matrix(rnorm(n_loc*n_var*(n_time+25)), n_loc*n_var))


y_true = aperm(array(tatato[-seq(n_loc*n_var*25)], c(n_var, n_loc, n_time)), c(2,1,3))
plot(rep(locs_no_na[,1], n_var), y_true[,,1], col = rep(seq(n_var), each = n_loc), cex = .3, pch = 15)
plot(rep(locs_no_na[,1], n_var), y_true[,,2], col = rep(seq(n_var), each = n_loc), cex = .3, pch = 15)
plot(rep(locs_no_na[,1], n_var), y_true[,,3], col = rep(seq(n_var), each = n_loc), cex = .3, pch = 15)
#plot(rep(locs_no_na[,1], n_var), y_true[,,100], col = rep(seq(n_var), each = n_loc), cex = .3, pch = 15)

# observed data set

y = y_true
##adding NAs
for(i in seq(50))y[1+floor(n_loc*runif(1)), 1+floor(n_var*runif(1)),] = NA
y[ceiling(runif(n = length(y)/100, max = length(y)))] = NA
y[5,,]=NA
#y[which((locs_no_na>.5)&(locs_no_na<1)),,50] = NA
#y[which((locs_no_na>1)&(locs_no_na<1.5)),,20] = NA

#adding noise
y = y + rnorm(length(y))


# dropping all NA idx
all_NA_idx = apply(y, c(1, 3), function(x)all(is.na(x)))
y = y[-which(all_NA_idx),,,drop = F]
locs = locs_no_na[-which(all_NA_idx),]

#plot(rep(locs[,1],3), y[,,50], col = rep(seq(3), each = nrow(locs)), pch  = 16, cex = .2)
#plot(rep(locs[,1],3), y[,,100], col = rep(seq(3), each = nrow(locs)), pch  = 16, cex = .2)

# markov chain states
Vecchia_approx_DAG = make_simple_Vecchia_approx_DAG(
  y = y, locs = locs)

vecchia_blocks = vecchia_block_approx(
  Vecchia_approx_DAG = Vecchia_approx_DAG, locs = locs, 
  time_depth = 3, #does not depend on params
  rho_vec = rho_vec, a = a, b = b, cc = cc, delta = delta, lambda = lambda, r = r, 
  A_vec = A_vec, nu_vec = nu_vec, log_range_vec = log_range_vec
)

X_noise = array(1, c(nrow(locs), 1, n_time))
X = array(1, c(nrow(locs), 1, n_time))
X_scale = array(1, c(nrow(locs), 1, n_time))

# initializing #########


mcmc_nngp_list = multivariate_NNGP_initialize(
  y = y, locs = locs, X = X, X_noise = X_noise, X_scale = X_scale, Vecchia_approx_DAG = Vecchia_approx_DAG, n_chains = 2, 
  time_depth = 5)
stuff_for_plots = list(
  locs_no_na = locs_no_na, 
  n_loc = n_loc,
  n_var = n_var, 
  y_true = y_true
)
#remove(list = setdiff(ls(), c("mcmc_nngp_list", "stuff_for_plots")))

source("MultiNNGP/R/GMA.R")
source("MultiNNGP/R/array_useful_functions.R")
source("MultiNNGP/R/multivariate_NN.R")
source("MultiNNGP/R/vecchia.R")
source("MultiNNGP/R/initialize.R")
source("MultiNNGP/R/latent_field_sampling.R")
source("MultiNNGP/R/basis_functions.R")
source("MultiNNGP/R/noise_variance.R")
#source("MultiNNGP/R/update.R")


chain = make_one_chain(mcmc_nngp_list)
chain$params$A_vec[]  = A_vec
chain$params$rho_vec[]  = rho_vec
chain$params$a_scal[] = a
chain$params$b_scal[] = b
chain$params$cc[]     = cc
chain$params$delta[]  = delta
chain$params$r[]      = r
chain$params$lambda[] = lambda



chain$params$log_range_vec[] = log_range_vec
chain$params$smoothness_vec[] = nu_vec

chain$stuff$vecchia_blocks = 
  vecchia_block_approx( 
    Vecchia_approx_DAG = mcmc_nngp_list$Vecchia_approx_DAG, locs = mcmc_nngp_list$locs, 
    time_depth = mcmc_nngp_list$useful_stuff$time_depth, #does not depend on params
    rho_vec = chain$params$rho_vec, a =  chain$params$a_scal,
    b = chain$params$b_scal, cc = chain$params$cc, delta = chain$params$delta, 
    lambda = chain$params$lambda, r = chain$params$r, 
    A_vec = chain$params$A_vec, nu_vec = chain$params$smoothness_vec, log_range_vec = chain$params$log_range_vec
  )

chain$stuff$transposed_vecchia_blocks = vecchia_blocks_t(chain$stuff$vecchia_blocks)
chain$stuff$precision_blocks = get_precision_blocks(vecchia_blocks = chain$stuff$vecchia_blocks)


# updating

#tatato = do_100_updates(chain, mcmc_nngp_list, kernel_learning_rate = 1, thinning_rate = 1)

