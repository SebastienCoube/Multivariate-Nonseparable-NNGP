##############
# Space only #
##############

source("MultiNNGP/R/GMA.R")
source("MultiNNGP/R/multivariate_NN.R")
source("MultiNNGP/R/vecchia.R")
source("MultiNNGP/R/initialize.R")
source("MultiNNGP/R/latent_field_sampling.R")

# space time layout
n_loc = 10000
n_var = 3
n_time = 1
locs_no_na = cbind(2*runif(n_loc), 1)


# parameters
rho_vec = rep(.8, n_var*(n_var-1)/2)
nu_vec = c(.5, 1.3, 2.5)
a2_vec = c(10,100,1000)


# sampling latent field
Vecchia_approx_DAG_no_NA = make_simple_Vecchia_approx_DAG(
  y = array(1, dim = c(n_loc, n_var, 1)), locs = locs_no_na
)
lower_tri_idx_DAG = get_lower_tri_idx_DAG(Vecchia_approx_DAG_no_NA$DAG)
var_idx = get_var_idx(Vecchia_approx_DAG_no_NA$DAG, Vecchia_approx_DAG_no_NA$field_position$var_idx)
vecchia_blocks_no_NA = vecchia_block_approx(
    Vecchia_approx_DAG = Vecchia_approx_DAG_no_NA, locs = locs_no_na, 
    lower_tri_idx_DAG = lower_tri_idx_DAG, 
    var_idx= var_idx, time_depth = 1, #does not depend on params
    rho_vec = rho_vec, a = NULL, b = NULL, cc = NULL, delta = NULL, lambda = NULL, r = NULL, 
    A_vec = NULL, nu_vec = nu_vec, a2_vec = a2_vec
)


y_true = array(dim = c(n_loc, n_var, 1))
y_true[] = t(matrix(vecchia_blocks_solve(vecchia_blocks = vecchia_blocks_no_NA, x = rnorm(n_loc*n_var)),3))

plot(rep(locs_no_na[,1],3), y_true, col = rep(seq(3), each = nrow(locs_no_na)), pch  = 16, cex = .5)

# observed data set
y = y_true + 1*rnorm(length(y_true))
y[ceiling(runif(n = 10000, max = length(y)))] = NA
all_NA_idx = apply(y, c(1, 3), function(x)all(is.na(x)))
y = y[-which(all_NA_idx),,,drop = F]
locs = locs_no_na[-which(all_NA_idx),]

plot(rep(locs[,1],3), y, col = rep(seq(3), each = nrow(locs)), pch  = 16, cex = .2)

# markov chain states
Vecchia_approx_DAG = make_simple_Vecchia_approx_DAG(
  y = y, locs = locs
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

X_noise = array(1, c(nrow(locs), 1, 1))
X = array(1, c(nrow(locs), 1, 1))
X_scale = array(1, c(nrow(locs), 1, 1))

mcmc_nngp_list = multivariate_NNGP_initialize(
  y = y, locs = locs, X = X, X_noise = X_noise, X_scale = X_scale, Vecchia_approx_DAG = Vecchia_approx_DAG, n_chains = 2)


precision_blocks = get_precision_blocks(vecchia_blocks)
noise_info = mcmc_nngp_list$chains$chain_1$stuff$noise_info
X_noise_list = mcmc_nngp_list$covariates$X_noise
diag_precision_block = precision_blocks[[1]]
useful_stuff = mcmc_nngp_list$useful_stuff


i=3
# subsetting locations
loc_subsetting = kmeans(locs, 5+rpois(1, 5))$clust
plot(locs[,1], loc_subsetting)
field_subsetting = loc_subsetting[mcmc_nngp_list$Vecchia_approx_DAG$field_position$location_idx]
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
  ##a_posteriori_diag_precision = get_a_posteriori_diag_precision(
  ##  diag_precision_block = diag_precision_block, 
  ##  useful_stuff = useful_stuff, 
  ##  loc_subset_idx = loc_subset_idx, 
  ##  noise_info = noise_info, 
  ##  X_noise_list = X_noise_list, 
  ##  Vecchia_approx_DAG = Vecchia_approx_DAG
  ##)
  field_obs_precision = get_field_obs_precision(
    noise_info = noise_info, 
    Vecchia_approx_DAG = Vecchia_approx_DAG, 
    loc_subset_idx = loc_subset_idx, 
    useful_stuff = useful_stuff, 
    X_noise_list = X_noise_list, 
    y = mcmc_nngp_list$y
  )
  
  mcmc_nngp_list$chains$chain_1$params$field[field_subset_idx] = 
    posterior_precision_solve_chol( # chol(QAA)-1
      a_posteriori_diag_precision_chols = a_posteriori_diag_precision_chols,
      below_diag_precision_blocks = NULL, # here NULL because 1 time period !
      vecchia_blocks = vecchia_blocks, 
      loc_subset_idx = loc_subset_idx, 
      x = rnorm(length(field_subset_idx)) - # normal sample following conditional variance
        # + conditional mean
        posterior_precision_solve_t_chol(# chol(QAA)-T -> combines with chol(QAA)-1 to give QAA -1
          a_posteriori_diag_precision_chols = a_posteriori_diag_precision_chols, 
          below_diag_precision_blocks = NULL, # here NULL because 1 time period !
          vecchia_blocks = vecchia_blocks, 
          loc_subset_idx = loc_subset_idx, 
          x = # QAB XB-muB
            (as.vector(precision_blocks[[1]]%*%(c(mcmc_nngp_list$chains$chain_1$params$field) * (field_subsetting != i))))[field_subset_idx] + 
            - as.vector(field_obs_precision[[1]] %*% c(t(mcmc_nngp_list$useful_stuff$y_na_killed[loc_subset_idx,,])))
        )
    )
  
  
  
  plot(rep(locs_no_na[,1],3), y_true, col = rep(seq(3), each = nrow(locs_no_na)), pch  = 16, cex = .3)
  points(mcmc_nngp_list$locs[Vecchia_approx_DAG$field_position$location_idx, 1], mcmc_nngp_list$chains$chain_1$params$field,  col = Vecchia_approx_DAG$field_position$var_idx, pch = 3, cex = .3)
  
}











##############
# Space time #
##############


source("MultiNNGP/R/GMA.R")
source("MultiNNGP/R/multivariate_NN.R")
source("MultiNNGP/R/vecchia.R")
source("MultiNNGP/R/initialize.R")
source("MultiNNGP/R/latent_field_sampling.R")



# space time layout
n_loc = 5000
n_var = 3
n_time = 150
locs_no_na = cbind(2*runif(n_loc), 1)


# parameters
rho_vec = rep(.8, n_var*(n_var-1)/2)
nu_vec = c(.5, 1.3, 2.5)
a2_vec = c(10,100,1000)
A_vec = rep(.5, n_var)
a =       .5
b =       1
delta=    .5
lambda =  .5

cc =      .05
r =       .05

# sampling latent field
Vecchia_approx_DAG_no_NA = make_simple_Vecchia_approx_DAG(
  y = array(1, dim = c(n_loc, n_var, n_time)), locs = locs_no_na
)
lower_tri_idx_DAG = get_lower_tri_idx_DAG(Vecchia_approx_DAG_no_NA$DAG)
var_idx = get_var_idx(Vecchia_approx_DAG_no_NA$DAG, Vecchia_approx_DAG_no_NA$field_position$var_idx)
vecchia_blocks_no_NA = vecchia_block_approx(
    Vecchia_approx_DAG = Vecchia_approx_DAG_no_NA, locs = locs_no_na, 
    lower_tri_idx_DAG = lower_tri_idx_DAG, 
    var_idx= var_idx, time_depth = 5, #does not depend on params
    rho_vec = rho_vec, a = a, b = b, cc = cc, delta = delta, lambda = lambda, r = r, 
    A_vec = A_vec, nu_vec = nu_vec, a2_vec = a2_vec
)

tatato = vecchia_blocks_solve(vecchia_blocks = vecchia_blocks_no_NA, x = rnorm(n_loc*n_var*(n_time+10)))


y_true = aperm(array(tatato[-seq(n_loc*n_var*10)], c(n_var, n_loc, n_time)), c(2,1,3))
plot(rep(locs_no_na[,1], n_var), y_true[,,1], col = rep(seq(3), each = n_loc), cex = .3, pch = 15)
plot(rep(locs_no_na[,1], n_var), y_true[,,25], col = rep(seq(3), each = n_loc), cex = .3, pch = 15)
plot(rep(locs_no_na[,1], n_var), y_true[,,50], col = rep(seq(3), each = n_loc), cex = .3, pch = 15)
plot(rep(locs_no_na[,1], n_var), y_true[,,150], col = rep(seq(3), each = n_loc), cex = .3, pch = 15)

# observed data set

y = y_true
#adding NAs
for(i in seq(1500))y[1+floor(n_loc*runif(1)), 1+floor(n_var*runif(1)),] = NA
y[ceiling(runif(n = 75000, max = length(y)))] = NA
y[which((locs_no_na>.5)&(locs_no_na<1)),,50] = NA
y[which((locs_no_na>1)&(locs_no_na<1.5)),,100] = NA
y[which((locs_no_na>1.5)&(locs_no_na<2)),,150] = NA

#adding noise
y = y + rnorm(length(y))

# dropping all NA idx
all_NA_idx = apply(y, c(1, 3), function(x)all(is.na(x)))
y = y[-which(all_NA_idx),,,drop = F]
locs = locs_no_na[-which(all_NA_idx),]

plot(rep(locs[,1],3), y[,,50], col = rep(seq(3), each = nrow(locs)), pch  = 16, cex = .2)
plot(rep(locs[,1],3), y[,,100], col = rep(seq(3), each = nrow(locs)), pch  = 16, cex = .2)
plot(rep(locs[,1],3), y[,,150], col = rep(seq(3), each = nrow(locs)), pch  = 16, cex = .2)

# markov chain states
Vecchia_approx_DAG = make_simple_Vecchia_approx_DAG(
  y = y, locs = locs
)
lower_tri_idx_DAG = get_lower_tri_idx_DAG(Vecchia_approx_DAG$DAG)
var_idx = get_var_idx(Vecchia_approx_DAG$DAG, Vecchia_approx_DAG$field_position$var_idx)
vecchia_blocks = vecchia_block_approx(
  Vecchia_approx_DAG = Vecchia_approx_DAG, locs = locs, 
  lower_tri_idx_DAG = lower_tri_idx_DAG, 
  var_idx= var_idx, time_depth = 5, #does not depend on params
  rho_vec = rho_vec, a = a, b = b, cc = cc, delta = delta, lambda = lambda, r = r, 
  A_vec = A_vec, nu_vec = nu_vec, a2_vec = a2_vec
)

X_noise = array(1, c(nrow(locs), 1, n_time))
X = array(1, c(nrow(locs), 1, n_time))
X_scale = array(1, c(nrow(locs), 1, n_time))

mcmc_nngp_list = multivariate_NNGP_initialize(
  y = y, locs = locs, X = X, X_noise = X_noise, X_scale = X_scale, Vecchia_approx_DAG = Vecchia_approx_DAG, n_chains = 2)


precision_blocks = get_precision_blocks(vecchia_blocks)
noise_info = mcmc_nngp_list$chains$chain_1$stuff$noise_info
X_noise_list = mcmc_nngp_list$covariates$X_noise
#diag_precision_block = precision_blocks[[1]]
useful_stuff = mcmc_nngp_list$useful_stuff
time_begin = 6
time_end = 175

noise_precisions = get_noise_precisions(
  noise_info = noise_info, 
  useful_stuff = useful_stuff, 
  Vecchia_approx_DAG = Vecchia_approx_DAG, 
  time_begin = useful_stuff$time_depth, time_end = useful_stuff$n_time_periods
  )


# gross partition
n_clust = rpois(1, 10)
clust = kmeans(locs, n_clust)$cl

#  selected_field = unlist(Vecchia_approx_DAG$field_position$loc_match[selected_locs])
#  selected_var = Vecchia_approx_DAG$field_position$var_idx[selected_field]

# cluster functions
  B = Matrix::sparseMatrix(
    i = selected_field, 
    j = selected_var, 
    dims = c(useful_stuff$n_field, useful_stuff$n_var_y)
  )
  B = Matrix::sparseMatrix(
    i = seq(useful_stuff$n_field), 
    j = useful_stuff$n_var_y * (clust[Vecchia_approx_DAG$field_position$location_idx[seq(useful_stuff$n_field)]]-1) + Vecchia_approx_DAG$field_position$var_idx, 
    dims = c(useful_stuff$n_field, max(clust)*useful_stuff$n_var_y)
  )
  
  
  # prior
  Q_B = lapply(precision_blocks, function(x)x %*% B)
  B_Q_B = lapply(Q_B, function(x)as.matrix(Matrix::crossprod(B,x)))
  prior_precision = get_block_toeplitz(B_Q_B, useful_stuff$n_time_periods- useful_stuff$time_depth+1)
  Q_field = vecchia_blocks_t_mult(
    x = vecchia_blocks_mult(
      x = mcmc_nngp_list$chains$chain_1$params$field,
      vecchia_blocks = vecchia_blocks
    ),
    vecchia_blocks = vecchia_blocks
  )[-seq(useful_stuff$n_field*(useful_stuff$time_depth-1))]
  B_Q_field = c(as.matrix(Matrix::crossprod(B, matrix(Q_field, useful_stuff$n_field))))
  
  # likelihood  
  tau_B = lapply(noise_precisions, function(x)x%*%B)
  B_tau_B = lapply(tau_B, function(x)as.matrix(Matrix::crossprod(B,x)))
  y_minus_field = mcmc_nngp_list$chains$chain_1$params$field[,useful_stuff$time_depth:useful_stuff$n_time_periods] - useful_stuff$y_loc_var_format_no_NA[,useful_stuff$time_depth:useful_stuff$n_time_periods]
  B_tau_y_minus_field = unlist(mapply(function(x,y) as.vector(Matrix::crossprod(x, y)), tau_B, split(y_minus_field, col(y_minus_field)), SIMPLIFY = F))
  
  # posterior
  posterior_precision_chol = 
    Matrix::expand(
      Matrix::Cholesky(
        Matrix::bdiag(B_tau_B) + 
          Matrix::sparseMatrix(
            i = prior_precision$i, 
            j = prior_precision$j, 
            x = prior_precision$x, 
            symmetric = T
          )
      )
    )
  
  # sampling
  eps = Matrix::t(posterior_precision_chol$P) %*%
    Matrix::solve(posterior_precision_chol$L, 
                  Matrix::solve(
                    Matrix::t(posterior_precision_chol$L), 
                    posterior_precision_chol$P %*% (-B_Q_field - B_tau_y_minus_field))
    )

  mcmc_nngp_list$chains$chain_1$params$field[,-seq(useful_stuff$time_depth-1)]  = mcmc_nngp_list$chains$chain_1$params$field[,-seq(useful_stuff$time_depth-1)] +
    as.matrix(B %*% matrix(eps, ncol = useful_stuff$n_time_periods - useful_stuff$time_depth+1))


plot(rep(locs_no_na[,1], n_var), y_true[,,1], col = rep(seq(3), each = n_loc), cex = .3, pch = 15)
points(mcmc_nngp_list$locs[Vecchia_approx_DAG$field_position$location_idx,1], 
       mcmc_nngp_list$chains$chain_1$params$field[,useful_stuff$buffer_depth], 
       cex = .3, pch = 16, col = Vecchia_approx_DAG$field_position$var_idx
       )

plot(rep(locs_no_na[,1], n_var), y_true[,,25], col = rep(seq(3), each = n_loc), cex = .3, pch = 15)
points(mcmc_nngp_list$locs[Vecchia_approx_DAG$field_position$location_idx,1], 
       mcmc_nngp_list$chains$chain_1$params$field[,useful_stuff$buffer_depth+25], 
       cex = .3, pch = 16, col = Vecchia_approx_DAG$field_position$var_idx
       )

plot(rep(locs_no_na[,1], n_var), y_true[,,50], col = rep(seq(3), each = n_loc), cex = .3, pch = 15)
points(mcmc_nngp_list$locs[Vecchia_approx_DAG$field_position$location_idx,1], 
       mcmc_nngp_list$chains$chain_1$params$field[,useful_stuff$buffer_depth+50], 
       cex = .3, pch = 16, col = Vecchia_approx_DAG$field_position$var_idx
       )

plot(rep(locs_no_na[,1], n_var), y_true[,,150], col = rep(seq(3), each = n_loc), cex = .3, pch = 15)
points(mcmc_nngp_list$locs[Vecchia_approx_DAG$field_position$location_idx,1], 
       mcmc_nngp_list$chains$chain_1$params$field[,useful_stuff$buffer_depth+150], 
       cex = .3, pch = 16, col = Vecchia_approx_DAG$field_position$var_idx
       )
  
  
  
  
  
n_impulsion_centers = 100
impulsion_centers = mcmc_nngp_list$locs[order(runif(useful_stuff$n_loc)<=n_impulsion_centers),]
impulsion_range = .2

impulsions = lapply(impulsion_centers, function(x)c(fields::rdist(x, mcmc_nngp_list$locs)))
impulsions = lapply(impulsions, function(x)
{
  i_idx = lapply(seq(1, 3), function(i)which(x<sqrt(5*impulsion_range/4^i)))
})
plot(locs[,1], rep(1, nrow(locs)))
points(impulsion_centers[1,1], 1, col=2)
points(locs[impulsions[[1]][[1]],1], rep(.8, length(impulsions[[1]][[1]])), col=2)
points(locs[impulsions[[1]][[2]],1], rep(.7, length(impulsions[[1]][[2]])), col=2)


# subsetting locations

time_window = c(70, 80)

for(i in seq(100))
{
  print(i)
  t1 = Sys.time()
  
  
  t1 = Sys.time()
  field_innovation = matrix(
    vecchia_blocks_solve(
      vecchia_blocks = vecchia_blocks, 
      x = rnorm((diff(time_window) + mcmc_nngp_list$useful_stuff$buffer_depth)*
                  mcmc_nngp_list$useful_stuff$n_field)),
    nrow = mcmc_nngp_list$useful_stuff$n_field
  )[,-seq(mcmc_nngp_list$useful_stuff$buffer_depth)]
  Sys.time()-t1
  
  t1 = Sys.time()
  
  #plot(mcmc_nngp_list$locs[mcmc_nngp_list$Vecchia_approx_DAG$field_position$location_idx,1],field_innovation[,1])
  #plot(mcmc_nngp_list$locs[mcmc_nngp_list$Vecchia_approx_DAG$field_position$location_idx,1],field_innovation[,10])
  #plot(locs[,1], space_tapering)
  field_innovation = field_innovation %*% Matrix::Diagonal(n = diff(time_window), x = time_tapering)
  #field_innovation = 
  #  Matrix::sparseMatrix(
  #  x = space_tapering[space_tapering>.01], 
  #  i = which(space_tapering>.01), 
  #  j = which(space_tapering>.01), 
  #  dims = c(length(space_tapering),length(space_tapering))
  #) %*% field_innovation
  time_tapering = exp(-((seq(diff(time_window)) -mean(seq(diff(time_window))))/(.25*diff(time_window)))^2)
  
  #plot(mcmc_nngp_list$locs[mcmc_nngp_list$Vecchia_approx_DAG$field_position$location_idx,1],field_innovation[,1], col = Vecchia_approx_DAG$field_position$var_idx)
  #plot(mcmc_nngp_list$locs[mcmc_nngp_list$Vecchia_approx_DAG$field_position$location_idx,1],field_innovation[,5], col = Vecchia_approx_DAG$field_position$var_idx)
  #plot(mcmc_nngp_list$locs[mcmc_nngp_list$Vecchia_approx_DAG$field_position$location_idx,1],field_innovation[,10], col = Vecchia_approx_DAG$field_position$var_idx)
  
  space_tapering = 
    exp(-.5*5000*fields::rdist(locs[1,,drop=F], locs)^6)
  selected_sites = which(space_tapering>.01)
  
  space_var_tapering = space_tapering[Vecchia_approx_DAG$field_position$location_idx]
  selected_space_var_pairs = which(space_var_tapering>.01)
  tapering_mat = Matrix::sparseMatrix(
    x = space_var_tapering[selected_space_var_pairs], 
    i = selected_space_var_pairs,
    j = rep(1, length(selected_space_var_pairs)),
    dims = c(nrow(field_innovation), 1)
  ) %*% as(t(time_tapering), "sparseMatrix")
Sys.time()-t1

t1 = Sys.time()
  proposed_field = 
    mcmc_nngp_list$chains$chain_1$params$field[,seq(time_window[1]-mcmc_nngp_list$useful_stuff$time_depth+1, 
                                                    time_window[2] + mcmc_nngp_list$useful_stuff$time_depth-1)]
  
  proposed_field[,seq(useful_stuff$time_depth, useful_stuff$time_depth+diff(time_window)-1)] = 
    as.matrix(
    proposed_field[,seq(useful_stuff$time_depth, useful_stuff$time_depth+ diff(time_window)-1)] * (1-.1*tapering_mat)^2 + 
    field_innovation * (.1*tapering_mat)^2
    )
Sys.time()-t1
  
  if(
    (
      -.5*sum(vecchia_blocks_mult(
        x = proposed_field,
        vecchia_blocks = vecchia_blocks
      )^2)
      +.5*sum(vecchia_blocks_mult(
        x = mcmc_nngp_list$chains$chain_1$params$field[,seq(time_window[1]-mcmc_nngp_list$useful_stuff$time_depth+1, 
                                                            time_window[2] + mcmc_nngp_list$useful_stuff$time_depth-1)],
        vecchia_blocks = vecchia_blocks
      )^2)
      + get_field_obs_density(
        selected_sites = selected_sites, 
        time_window = time_window, 
        field_subset = proposed_field[,seq(useful_stuff$time_depth, useful_stuff$time_depth + diff(time_window))], 
        y_ = mcmc_nngp_list$y, noise_info = noise_info, useful_stuff = useful_stuff, 
        Vecchia_approx_DAG = Vecchia_approx_DAG
      )       - get_field_obs_density(
        selected_sites = selected_sites, 
        time_window = time_window, 
        field_subset = mcmc_nngp_list$chains$chain_1$params$field[,seq(time_window[1], time_window[2])], 
        y_ = mcmc_nngp_list$y, noise_info = noise_info, useful_stuff = useful_stuff, 
        Vecchia_approx_DAG = Vecchia_approx_DAG
      )
    )
    > log(runif(1))
  )
  {
    mcmc_nngp_list$chains$chain_1$params$field[,seq(time_window[1], time_window[2])] =
      proposed_field[,seq(useful_stuff$time_depth, useful_stuff$time_depth + diff(time_window))]
    print("tatatoooo!")
  }
  
  Sys.time()-t1
  
  
  
  
}
plot(rep(locs_no_na[,1],3), y_true[,,50],  col = rep(seq(3), each = nrow(locs_no_na)), pch  = 16, cex = .2)
points(mcmc_nngp_list$locs[Vecchia_approx_DAG$field_position$location_idx, 1], mcmc_nngp_list$chains$chain_1$params$field[,75],  col = Vecchia_approx_DAG$field_position$var_idx, pch = 3, cex = .3)
















