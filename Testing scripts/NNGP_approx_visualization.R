source("MultiNNGP/R/GMA.R")
source("MultiNNGP/R/multivariate_NN.R")
# generate test data
set.seed(10)
n_var = 3

locs_1 = cbind(rep(seq(100)/100, 3), 1)
var_tag_1 = rep(seq(3), each= 100)

rho_vec =c(.8, .8, .8)

a2_vec = c(100, 200, 400)
nu_vec = c(1.5, 1, .5)
alpha = .5
a =  .9
b =  1
cc = .1
lambda = .5
delta = .5
r = .5
A_vec = rep(.5, 3)
u = seq(0, 50)


multiplier = get_multiplier(
  a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
  r = r, A_vec = A_vec, nu_vec = nu_vec, a2_vec = a2_vec, u = u
)
effective_range = get_effective_range(
  a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
  r = r, A_vec = A_vec, a2_vec = a2_vec, u = u
)

#### # carry on with toy example
covmat_coeffs_1  =   GMA_compressed(
  locs_ = locs_1, 
  lower_tri_idx = get_lower_tri_idx(nrow(locs_1), diag = T),
  var_idx = position_in_lower_tri_cross_vec(var_tag_1, n_var, diag = T),
  multiplier = multiplier#[1,,drop=FALSE]
  , 
  effective_range = effective_range,
  nu_ = expand_nu(nu_vec),
  rho_vec_with_ones = put_ones_in_rho_vec(rho_vec)
)


covmats = expand_covmat_into_blocks(covmat_coeffs = covmat_coeffs_1$covmat, n_loc = nrow(locs_1), block_lower_tri_idx = get_lower_tri_idx(nrow(locs_1), T))

plot(locs_1[,1], covmats[[1]][,50], ylim = c(0, 1), col = rep(seq(3), each = 100), pch = 16, cex = .5)
plot(locs_1[,1], covmats[[1]][,150], ylim = c(0, 1), col = rep(seq(3), each = 100), pch = 16, cex = .5)
plot(locs_1[,1], covmats[[1]][,250], ylim = c(0, 1), col = rep(seq(3), each = 100), pch = 16, cex = .5)



plot(locs_1[seq(100),1], covmats[[1]][seq(100),50], ylim = c(0, 1), pch = 16, cex = .5)
points(locs_1[seq(100),1], covmats[[5]][seq(100),50], ylim = c(0, 1), pch = 16, cex = .5)
points(locs_1[seq(100),1], covmats[[10]][seq(100),50], ylim = c(0, 1), pch = 16, cex = .5)


y = array(1, c(100, 3, 100))
Vecchia_approx_DAG = make_simple_Vecchia_approx_DAG(
  y = y, locs = locs_1
)


vecchia_blocks = vecchia_block_approx(
  Vecchia_approx_DAG = Vecchia_approx_DAG, locs = locs_1, 
  lower_tri_idx_DAG = lower_tri_idx_DAG, 
  var_idx= var_idx, time_depth = 5, #does not depend on params
  rho_vec = rho_vec, a = a, b = b, cc = cc, delta = delta, lambda = lambda, r = r, 
  A_vec = A_vec, nu_vec = nu_vec, log_range_vec = log_range_vec
)

