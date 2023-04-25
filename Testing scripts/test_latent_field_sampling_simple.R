##############
# Space time #
##############


source("MultiNNGP/R/GMA.R")
source("MultiNNGP/R/multivariate_NN.R")
source("MultiNNGP/R/vecchia.R")
source("MultiNNGP/R/initialize.R")
source("MultiNNGP/R/latent_field_sampling.R")
source("MultiNNGP/R/basis_functions.R")


# synthetic data set #########

# space time layout
n_loc = 3000
n_var = 3
n_time = 100
locs_no_na = cbind(2*runif(n_loc), 1)


# parameters
rho_vec = rep(.8, n_var*(n_var-1)/2)
nu_vec = c(.5, 1, 1.5)
a2_vec = c(10,30,50)
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
t1 = Sys.time()
vecchia_blocks_no_NA = vecchia_block_approx(
  Vecchia_approx_DAG = Vecchia_approx_DAG_no_NA, locs = locs_no_na, 
  lower_tri_idx_DAG = lower_tri_idx_DAG, 
  var_idx= var_idx, time_depth = 5, #does not depend on params
  rho_vec = rho_vec, a = a, b = b, cc = cc, delta = delta, lambda = lambda, r = r, 
  A_vec = A_vec, nu_vec = nu_vec, a2_vec = a2_vec
)
Sys.time()-t1

tatato = vecchia_blocks_solve(vecchia_blocks = vecchia_blocks_no_NA, 
                              x = matrix(rnorm(n_loc*n_var*(n_time+10)), n_loc*n_var))


y_true = aperm(array(tatato[-seq(n_loc*n_var*10)], c(n_var, n_loc, n_time)), c(2,1,3))
plot(rep(locs_no_na[,1], n_var), y_true[,,1], col = rep(seq(3), each = n_loc), cex = .3, pch = 15)
plot(rep(locs_no_na[,1], n_var), y_true[,,25], col = rep(seq(3), each = n_loc), cex = .3, pch = 15)
plot(rep(locs_no_na[,1], n_var), y_true[,,50], col = rep(seq(3), each = n_loc), cex = .3, pch = 15)

# observed data set

y = y_true
#adding NAs
for(i in seq(1500))y[1+floor(n_loc*runif(1)), 1+floor(n_var*runif(1)),] = NA
y[ceiling(runif(n = 75000, max = length(y)))] = NA
y[which((locs_no_na>.5)&(locs_no_na<1)),,50] = NA
y[which((locs_no_na>1)&(locs_no_na<1.5)),,100] = NA

#adding noise
y = y + rnorm(length(y))

# dropping all NA idx
all_NA_idx = apply(y, c(1, 3), function(x)all(is.na(x)))
y = y[-which(all_NA_idx),,,drop = F]
locs = locs_no_na[-which(all_NA_idx),]

plot(rep(locs[,1],3), y[,,50], col = rep(seq(3), each = nrow(locs)), pch  = 16, cex = .2)
plot(rep(locs[,1],3), y[,,100], col = rep(seq(3), each = nrow(locs)), pch  = 16, cex = .2)

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
transposed_vecchia_blocks = vecchia_blocks_t(vecchia_blocks = vecchia_blocks)

X_noise = array(1, c(nrow(locs), 1, n_time))
X = array(1, c(nrow(locs), 1, n_time))
X_scale = array(1, c(nrow(locs), 1, n_time))

# initializing #########


mcmc_nngp_list = multivariate_NNGP_initialize(
  y = y, locs = locs, X = X, X_noise = X_noise, X_scale = X_scale, Vecchia_approx_DAG = Vecchia_approx_DAG, n_chains = 2, 
  time_depth = 5)


precision_blocks = get_precision_blocks(vecchia_blocks)
noise_info = mcmc_nngp_list$chains$chain_1$stuff$noise_info
X_noise_list = mcmc_nngp_list$covariates$X_noise
#diag_precision_block = precision_blocks[[1]]
useful_stuff = mcmc_nngp_list$useful_stuff

time_begin = useful_stuff$time_depth
time_end = useful_stuff$n_time_periods - useful_stuff$time_depth +1

t1 = Sys.time()
noise_precisions = get_noise_precisions(
  noise_info = noise_info, 
  useful_stuff = useful_stuff, 
  Vecchia_approx_DAG = Vecchia_approx_DAG, 
  time_begin = time_begin, time_end = time_end
)
Sys.time()-t1

# using blocks 'n' bases #####

mcmc_nngp_list$chains$chain_1$params$field = 0*mcmc_nngp_list$chains$chain_1$params$field
record = NULL

##### 

image_idx = 0 
#pdf(paste("Gifs_field_sampling/blocks_and_bases/bnb_0.pdf"))
plot(rep(locs_no_na[,1], n_var), y_true[,,10], col = rep(c("lightgray", "lightpink", "lightblue"), each = n_loc), cex = .3, pch = 15, 
     xlab = "spatial site", 
     ylab = "latent field and true field", 
     main = "initial state"
)
legend("topright", 
       legend = c(
         "true field 1"   , 
         "true field 2"   , 
         "true field 3"   , 
         "sampled field 1",
         "sampled field 2",
         "sampled field 3"
         ), 
       fill = c(
          "lightgray",
          "lightpink"  ,
          "lightblue" ,
          "black",
          "red"  ,
          "blue" 
       )
       )
points(mcmc_nngp_list$locs[Vecchia_approx_DAG$field_position$location_idx,1], 
       mcmc_nngp_list$chains$chain_1$params$field[,useful_stuff$buffer_depth + 10], 
       cex = .3, pch = 1, col = Vecchia_approx_DAG$field_position$var_idx + (Vecchia_approx_DAG$field_position$var_idx==3)
)
#dev.off()



for(i in seq(5))
{
  cluster_size_target = 100
  locs_ = locs[Vecchia_approx_DAG$field_position$location_idx,]
  
  recursive_k_means_clust = recursive_k_means(locs_ = locs[Vecchia_approx_DAG$field_position$location_idx,], cluster_size_target)
  hist(table(recursive_k_means_clust$clust))
  
  basis_functions = c(
    get_indicator_basis_functions(recursive_k_means_clust = recursive_k_means_clust, useful_stuff = useful_stuff)
    , 
    unlist(lapply(
      get_grids(
        points = locs_, 
        cluster_size_target = cluster_size_target), 
      function(x)get_basis_functions(points = locs_, tile_info = x, cluster_size_target = cluster_size_target, Vecchia_approx_DAG = Vecchia_approx_DAG, useful_stuff = useful_stuff)), recursive = F)
  )
  
  coloring=  color_basis_functions(basis_functions)
  table(coloring)
  length(basis_functions)
  
  #plot(locs_[,1], basis_functions[[length(basis_functions)]][,1])
  #plot(locs_[,1], basis_functions[[length(basis_functions)-1]][,1])
  #plot(locs_[,1], basis_functions[[length(basis_functions)-2]][,1])
  #plot(locs_[,1], basis_functions[[length(basis_functions)-3]][,1])
  #plot(locs_[,1], basis_functions[[length(basis_functions)-4]][,4])
  
  #Something not working in the solving iot sample the field
  
  
  # sampling #####
  colo = 3
  
  
  ddddd = Sys.time()
  for(colo in seq(max(coloring)))
  { 
    print(colo)
    t1 = Sys.time()
    # precision times latent field
    Q_field = vecchia_blocks_mult(
      x = mcmc_nngp_list$chains$chain_1$params$field,
      vecchia_blocks = vecchia_blocks
    )
    Q_field[,-seq(time_begin, time_end)] = 0
    Q_field = vecchia_blocks_t_mult(
      x = Q_field,
      transposed_vecchia_blocks = transposed_vecchia_blocks
    )
    Q_field = matrix(Q_field, nrow = useful_stuff$n_field)
    # tau y - field
    y_minus_field = mcmc_nngp_list$chains$chain_1$params$field[,seq(time_begin, time_end)] - useful_stuff$y_loc_var_format_no_NA[,seq(time_begin, time_end)]
    tau_y_minus_field = do.call(cbind, mapply(function(x,y) as.vector(Matrix::crossprod(x, y)), noise_precisions, split(y_minus_field, col(y_minus_field)), SIMPLIFY = F))
    Sys.time()-t1
    
    t1 = Sys.time()
    i_basis_function=1
    eps = parallel::mclapply(
      mc.cores = 5, which(coloring == colo)
      ,
      function(i_basis_function)
      {
        #      # part with observations ####
        B_tau_y_minus_field = Matrix::crossprod(basis_functions[[i_basis_function]], tau_y_minus_field)
        
        B_tau_B_x = list()
        B_tau_B_i = list()
        B_tau_B_j = list()
        for(i_time in seq(time_begin, time_end)){
          B_tau_B = Matrix::crossprod(basis_functions[[i_basis_function]], Matrix::crossprod(noise_precisions[[i_time-time_begin+1]], basis_functions[[i_basis_function]]))
          B_tau_B = Matrix::forceSymmetric(B_tau_B)
          if(length(B_tau_B@x)>0)
          {
            B_tau_B_x = c(B_tau_B_x, list(B_tau_B@x))
            B_tau_B_i = c(B_tau_B_i, list((i_time-time_begin)*ncol(basis_functions[[i_basis_function]]) + B_tau_B@i+1))
            B_tau_B_j = c(B_tau_B_j, list((i_time-time_begin)*ncol(basis_functions[[i_basis_function]]) + findInterval(seq(length(B_tau_B@x))-1,B_tau_B@p[-1])+1))
          }
        }
        B_tau_B_x = unlist(B_tau_B_x)
        B_tau_B_i = unlist(B_tau_B_i)
        B_tau_B_j = unlist(B_tau_B_j)
        
        #Sys.time()-t1
        # part with prior ####
        B_Q_B = lapply(precision_blocks, function(x)Matrix::crossprod(basis_functions[[i_basis_function]], x) %*% basis_functions[[i_basis_function]])
        B_Q_B[[1]]=Matrix::forceSymmetric(B_Q_B[[1]], uplo = "L")
        prior_precision = get_block_toeplitz(B_Q_B, time_end- time_begin+1)
        B_Q_field = (Matrix::crossprod(basis_functions[[i_basis_function]], Q_field[,seq(time_begin, time_end)]))
        
        # posterior ####
        t1 = Sys.time()
        posterior_chol = Matrix::chol(Matrix::sparseMatrix(
          i = c(prior_precision$i, B_tau_B_j), 
          j = c(prior_precision$j, B_tau_B_i), 
          x = c(prior_precision$x, B_tau_B_x), 
          symmetric = T, 
          dims = rep((time_end-time_begin+1)*ncol(basis_functions[[i_basis_function]]), 2)
        ))
        Sys.time()-t1
        
        # sampling coefficients ####
        
        eps = 
          matrix(
            Matrix::solve(
              posterior_chol, 
              rnorm(length(B_tau_y_minus_field))
              + Matrix::solve(
                Matrix::t(posterior_chol), 
                as.vector(
                  - B_Q_field 
                  - B_tau_y_minus_field
                )
              )
            ), 
            nrow = ncol(basis_functions[[i_basis_function]])
          )
        
        return(eps)
      }, mc.preschedule = F
    )
    Sys.time()-t1
    
    mcmc_nngp_list$chains$chain_1$params$field[,time_begin:time_end] = 
      mcmc_nngp_list$chains$chain_1$params$field[,time_begin:time_end] +
      as.matrix(Reduce("+", mapply(function(e, b) b%*%e, eps, basis_functions[which(coloring==colo)])))
    
    
    
    
    
    image_idx = image_idx + 1
    #pdf(paste("Gifs_field_sampling/blocks_and_bases/bnb_", sep = "", image_idx, ".pdf"))
    plot(rep(locs_no_na[,1], n_var), y_true[,,10], col = rep(c("lightgray", "lightpink", "lightblue"), each = n_loc), cex = .3, pch = 15, 
         xlab = "spatial site", 
         ylab = "latent field and true field", 
         main = paste("pass", i, "color", colo)
    )
    legend("topright", 
           legend = c(
             "true field 1"   , 
             "true field 2"   , 
             "true field 3"   , 
             "sampled field 1",
             "sampled field 2",
             "sampled field 3"
           ), 
           fill = c(
             "lightgray",
                        "lightpink"  ,
                        "lightblue" ,
                        "black",
                        "red"  ,
                        "blue"
           )
    )
    points(mcmc_nngp_list$locs[Vecchia_approx_DAG$field_position$location_idx,1], 
           mcmc_nngp_list$chains$chain_1$params$field[,useful_stuff$buffer_depth + 10], 
           cex = .3, pch = 1, col = Vecchia_approx_DAG$field_position$var_idx + (Vecchia_approx_DAG$field_position$var_idx==3)
    )
    #dev.off()
    
    
    
  }
  record = cbind(record, mcmc_nngp_list$chains$chain_1$params$field[,useful_stuff$buffer_depth + 10])
  fffff = Sys.time()-ddddd
}


# using only blocks #####


mcmc_nngp_list$chains$chain_1$params$field = 0*mcmc_nngp_list$chains$chain_1$params$field

image_idx = 0

pdf(paste("Gifs_field_sampling/only_blocks/ob_0.pdf"))
plot(rep(locs_no_na[,1], n_var), y_true[,,10], col = rep(c("lightgray", "lightpink", "lightblue"), each = n_loc), cex = .3, pch = 15, 
     xlab = "spatial site", 
     ylab = "latent field and true field", 
     main = "initial state"
)
legend("topright", 
       legend = c(
         "true field 1"   , 
         "true field 2"   , 
         "true field 3"   , 
         "sampled field 1",
         "sampled field 2",
         "sampled field 3"
       ), 
       fill = c(
         "lightgray",
                    "lightpink"  ,
                    "lightblue" ,
                    "black",
                    "red"  ,
                    "blue" 
       )
)
points(mcmc_nngp_list$locs[Vecchia_approx_DAG$field_position$location_idx,1], 
       mcmc_nngp_list$chains$chain_1$params$field[,useful_stuff$buffer_depth + 10], 
       cex = .3, pch = 1, col = Vecchia_approx_DAG$field_position$var_idx + (Vecchia_approx_DAG$field_position$var_idx==3)
)
dev.off()



for(i in seq(10))
{
  cluster_size_target = 100
  locs_ = locs[Vecchia_approx_DAG$field_position$location_idx,]
  
  recursive_k_means_clust = recursive_k_means(locs_ = locs[Vecchia_approx_DAG$field_position$location_idx,], cluster_size_target)
  hist(table(recursive_k_means_clust$clust))
  
  basis_functions = c(
    get_indicator_basis_functions(recursive_k_means_clust = recursive_k_means_clust, useful_stuff = useful_stuff)
  )
  
  coloring=  color_basis_functions(basis_functions)
  table(coloring)
  length(basis_functions)
  
  #plot(locs_[,1], basis_functions[[length(basis_functions)]][,1])
  #plot(locs_[,1], basis_functions[[length(basis_functions)-1]][,1])
  #plot(locs_[,1], basis_functions[[length(basis_functions)-2]][,1])
  #plot(locs_[,1], basis_functions[[length(basis_functions)-3]][,1])
  #plot(locs_[,1], basis_functions[[length(basis_functions)-4]][,4])
  
  #Something not working in the solving iot sample the field
  
  
  # sampling #####
  colo = 3
  
  
  ddddd = Sys.time()
  for(colo in seq(max(coloring)))
  { 
    print(colo)
    t1 = Sys.time()
    # precision times latent field
    Q_field = vecchia_blocks_mult(
      x = mcmc_nngp_list$chains$chain_1$params$field,
      vecchia_blocks = vecchia_blocks
    )
    Q_field[,-seq(time_begin, time_end)] = 0
    Q_field = vecchia_blocks_t_mult(
      x = Q_field,
      transposed_vecchia_blocks = transposed_vecchia_blocks
    )
    Q_field = matrix(Q_field, nrow = useful_stuff$n_field)
    # tau y - field
    y_minus_field = mcmc_nngp_list$chains$chain_1$params$field[,seq(time_begin, time_end)] - useful_stuff$y_loc_var_format_no_NA[,seq(time_begin, time_end)]
    tau_y_minus_field = do.call(cbind, mapply(function(x,y) as.vector(Matrix::crossprod(x, y)), noise_precisions, split(y_minus_field, col(y_minus_field)), SIMPLIFY = F))
    Sys.time()-t1
    
    t1 = Sys.time()
    i_basis_function=1
    eps = parallel::mclapply(
      mc.cores = 5, which(coloring == colo)
      ,
      function(i_basis_function)
      {
        #      # part with observations ####
        B_tau_y_minus_field = Matrix::crossprod(basis_functions[[i_basis_function]], tau_y_minus_field)
        
        B_tau_B_x = list()
        B_tau_B_i = list()
        B_tau_B_j = list()
        for(i_time in seq(time_begin, time_end)){
          B_tau_B = Matrix::crossprod(basis_functions[[i_basis_function]], Matrix::crossprod(noise_precisions[[i_time-time_begin+1]], basis_functions[[i_basis_function]]))
          B_tau_B = Matrix::forceSymmetric(B_tau_B)
          if(length(B_tau_B@x)>0)
          {
            B_tau_B_x = c(B_tau_B_x, list(B_tau_B@x))
            B_tau_B_i = c(B_tau_B_i, list((i_time-time_begin)*ncol(basis_functions[[i_basis_function]]) + B_tau_B@i+1))
            B_tau_B_j = c(B_tau_B_j, list((i_time-time_begin)*ncol(basis_functions[[i_basis_function]]) + findInterval(seq(length(B_tau_B@x))-1,B_tau_B@p[-1])+1))
          }
        }
        B_tau_B_x = unlist(B_tau_B_x)
        B_tau_B_i = unlist(B_tau_B_i)
        B_tau_B_j = unlist(B_tau_B_j)
        
        #Sys.time()-t1
        # part with prior ####
        B_Q_B = lapply(precision_blocks, function(x)Matrix::crossprod(basis_functions[[i_basis_function]], x) %*% basis_functions[[i_basis_function]])
        B_Q_B[[1]]=Matrix::forceSymmetric(B_Q_B[[1]], uplo = "L")
        prior_precision = get_block_toeplitz(B_Q_B, time_end- time_begin+1)
        B_Q_field = (Matrix::crossprod(basis_functions[[i_basis_function]], Q_field[,seq(time_begin, time_end)]))
        
        # posterior ####
        t1 = Sys.time()
        posterior_chol = Matrix::chol(Matrix::sparseMatrix(
          i = c(prior_precision$i, B_tau_B_j), 
          j = c(prior_precision$j, B_tau_B_i), 
          x = c(prior_precision$x, B_tau_B_x), 
          symmetric = T, 
          dims = rep((time_end-time_begin+1)*ncol(basis_functions[[i_basis_function]]), 2)
        ))
        Sys.time()-t1
        
        # sampling coefficients ####
        
        eps = 
          matrix(
            Matrix::solve(
              posterior_chol, 
              rnorm(length(B_tau_y_minus_field))
              + Matrix::solve(
                Matrix::t(posterior_chol), 
                as.vector(
                  - B_Q_field 
                  - B_tau_y_minus_field
                )
              )
            ), 
            nrow = ncol(basis_functions[[i_basis_function]])
          )
        
        return(eps)
      }, mc.preschedule = F
    )
    Sys.time()-t1
    
    mcmc_nngp_list$chains$chain_1$params$field[,time_begin:time_end] = 
      mcmc_nngp_list$chains$chain_1$params$field[,time_begin:time_end] +
      as.matrix(Reduce("+", mapply(function(e, b) b%*%e, eps, basis_functions[which(coloring==colo)])))
    
    
    
    
    
    image_idx = image_idx+1
    pdf(paste("Gifs_field_sampling/only_blocks/ob_", sep = "", image_idx, ".pdf"))
    plot(rep(locs_no_na[,1], n_var), y_true[,,10], col = rep(c("lightgray", "lightpink", "lightblue"), each = n_loc), cex = .3, pch = 15, 
         xlab = "spatial site", 
         ylab = "latent field and true field", 
         main = paste("pass", i, "color", colo)
    )
    legend("topright", 
           legend = c(
             "true field 1"   , 
             "true field 2"   , 
             "true field 3"   , 
             "sampled field 1",
             "sampled field 2",
             "sampled field 3"
           ), 
           fill = c(
             "lightgray",
                        "lightpink"  ,
                        "lightblue" ,
                        "black",
                        "red"  ,
                        "blue"
           )
    )
    points(mcmc_nngp_list$locs[Vecchia_approx_DAG$field_position$location_idx,1], 
           mcmc_nngp_list$chains$chain_1$params$field[,useful_stuff$buffer_depth + 10], 
           cex = .3, pch = 1, col = Vecchia_approx_DAG$field_position$var_idx + (Vecchia_approx_DAG$field_position$var_idx==3)
    )
    dev.off()
    
    
  }
  fffff = Sys.time()-ddddd
}


##############
# Space time, with two space dims  #
##############


source("MultiNNGP/R/GMA.R")
source("MultiNNGP/R/multivariate_NN.R")
source("MultiNNGP/R/vecchia.R")
source("MultiNNGP/R/initialize.R")
source("MultiNNGP/R/latent_field_sampling.R")
source("MultiNNGP/R/basis_functions.R")


# synthetic data set #########

# space time layout
n_loc = 2000
n_var = 3
n_time = 100
locs_no_na = cbind(2*runif(n_loc), 2*runif(n_loc))


# parameters
rho_vec = rep(.8, n_var*(n_var-1)/2)
nu_vec = c(.5, 1, 1.5)
a2_vec = c(10,30,50)
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
t1 = Sys.time()
vecchia_blocks_no_NA = vecchia_block_approx(
  Vecchia_approx_DAG = Vecchia_approx_DAG_no_NA, locs = locs_no_na, 
  lower_tri_idx_DAG = lower_tri_idx_DAG, 
  var_idx= var_idx, time_depth = 5, #does not depend on params
  rho_vec = rho_vec, a = a, b = b, cc = cc, delta = delta, lambda = lambda, r = r, 
  A_vec = A_vec, nu_vec = nu_vec, a2_vec = a2_vec
)
Sys.time()-t1

tatato = vecchia_blocks_solve(vecchia_blocks = vecchia_blocks_no_NA, 
                              x = matrix(rnorm(n_loc*n_var*(n_time+10)), n_loc*n_var))


y_true = aperm(array(tatato[-seq(n_loc*n_var*10)], c(n_var, n_loc, n_time)), c(2,1,3))
plot(rep(locs_no_na[,1], n_var), y_true[,,1], col = rep(seq(3), each = n_loc), cex = .3, pch = 15)
plot(rep(locs_no_na[,1], n_var), y_true[,,25], col = rep(seq(3), each = n_loc), cex = .3, pch = 15)
plot(rep(locs_no_na[,1], n_var), y_true[,,50], col = rep(seq(3), each = n_loc), cex = .3, pch = 15)

# observed data set

y = y_true
#adding NAs
for(i in seq(1500))y[1+floor(n_loc*runif(1)), 1+floor(n_var*runif(1)),] = NA

#adding noise
y = y + rnorm(length(y))

# dropping all NA idx
all_NA_idx = apply(y, c(1, 3), function(x)all(is.na(x)))
y = y[-which(all_NA_idx),,,drop = F]
locs = locs_no_na[-which(all_NA_idx),]

plot(rep(locs[,1],3), y[,,50], col = rep(seq(3), each = nrow(locs)), pch  = 16, cex = .2)
plot(rep(locs[,1],3), y[,,100], col = rep(seq(3), each = nrow(locs)), pch  = 16, cex = .2)

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
transposed_vecchia_blocks = vecchia_blocks_t(vecchia_blocks = vecchia_blocks)

X_noise = array(1, c(nrow(locs), 1, n_time))
X = array(1, c(nrow(locs), 1, n_time))
X_scale = array(1, c(nrow(locs), 1, n_time))

# initializing #########


mcmc_nngp_list = multivariate_NNGP_initialize(
  y = y, locs = locs, X = X, X_noise = X_noise, X_scale = X_scale, Vecchia_approx_DAG = Vecchia_approx_DAG, n_chains = 2, 
  time_depth = 5)


precision_blocks = get_precision_blocks(vecchia_blocks)
noise_info = mcmc_nngp_list$chains$chain_1$stuff$noise_info
X_noise_list = mcmc_nngp_list$covariates$X_noise
#diag_precision_block = precision_blocks[[1]]
useful_stuff = mcmc_nngp_list$useful_stuff

time_begin = useful_stuff$time_depth
time_end = useful_stuff$n_time_periods - useful_stuff$time_depth +1

t1 = Sys.time()
noise_precisions = get_noise_precisions(
  noise_info = noise_info, 
  useful_stuff = useful_stuff, 
  Vecchia_approx_DAG = Vecchia_approx_DAG, 
  time_begin = time_begin, time_end = time_end
)
Sys.time()-t1

#clustering and basis functions #####

cluster_size_target = 100
locs_ = locs[Vecchia_approx_DAG$field_position$location_idx,]

recursive_k_means_clust = recursive_k_means(locs_ = locs[Vecchia_approx_DAG$field_position$location_idx,], cluster_size_target)
hist(table(recursive_k_means_clust$clust))

basis_functions = c(
  get_indicator_basis_functions(recursive_k_means_clust = recursive_k_means_clust, useful_stuff = useful_stuff)
  , 
  unlist(lapply(
    get_grids(
      points = locs_, 
      cluster_size_target = cluster_size_target), 
    function(x)get_basis_functions(points = locs_, tile_info = x, cluster_size_target = cluster_size_target, Vecchia_approx_DAG = Vecchia_approx_DAG, useful_stuff = useful_stuff)), recursive = F)
)



coloring=  color_basis_functions(basis_functions)
table(coloring)
length(basis_functions)

plot(locs_[,1], basis_functions[[length(basis_functions)]][,1])
plot(locs_[,1], basis_functions[[length(basis_functions)-1]][,1])
plot(locs_[,1], basis_functions[[length(basis_functions)-3]][,1])
plot(locs_[,1], basis_functions[[length(basis_functions)-5]][,1])
plot(locs_[,1], basis_functions[[length(basis_functions)-6]][,1])


# sampling #####
colo = 3
#mcmc_nngp_list$chains$chain_1$params$field = 0*mcmc_nngp_list$chains$chain_1$params$field



for(colo in seq(max(coloring)))
{ 
  
  plot(rep(locs_no_na[,1], n_var), y_true[,,50], col = rep(seq(3), each = n_loc), cex = .3, pch = 15)
  points(mcmc_nngp_list$locs[Vecchia_approx_DAG$field_position$location_idx,1], 
         mcmc_nngp_list$chains$chain_1$params$field[,useful_stuff$buffer_depth + 50], 
         cex = .3, pch = 1, col = Vecchia_approx_DAG$field_position$var_idx
  )
  
  plot(rep(locs_no_na[,1], n_var), y_true[,,1], col = rep(seq(3), each = n_loc), cex = .3, pch = 15)
  points(mcmc_nngp_list$locs[Vecchia_approx_DAG$field_position$location_idx,1], 
         mcmc_nngp_list$chains$chain_1$params$field[,useful_stuff$buffer_depth], 
         cex = .3, pch = 1, col = Vecchia_approx_DAG$field_position$var_idx
  )
  plot(rep(locs_no_na[,1], n_var), y_true[,,10], col = rep(c("lightgray", "lightpink", "lightblue"), each = n_loc), cex = .3, pch = 15, )
  points(mcmc_nngp_list$locs[Vecchia_approx_DAG$field_position$location_idx,1], 
         mcmc_nngp_list$chains$chain_1$params$field[,useful_stuff$buffer_depth + 10], 
         cex = .3, pch = 1, col = Vecchia_approx_DAG$field_position$var_idx + (Vecchia_approx_DAG$field_position$var_idx==3)
  )
  
  
  
  
  
  print(colo)
  
  t1 = Sys.time()
  # precision times latent field
  Q_field = vecchia_blocks_mult(
    x = mcmc_nngp_list$chains$chain_1$params$field,
    vecchia_blocks = vecchia_blocks
  )
  Q_field[,-seq(time_begin, time_end)] = 0
  Q_field = vecchia_blocks_t_mult(
    x = Q_field,
    transposed_vecchia_blocks = transposed_vecchia_blocks
  )
  Q_field = matrix(Q_field, nrow = useful_stuff$n_field)
  # tau y - field
  y_minus_field = mcmc_nngp_list$chains$chain_1$params$field[,seq(time_begin, time_end)] - useful_stuff$y_loc_var_format_no_NA[,seq(time_begin, time_end)]
  tau_y_minus_field = do.call(cbind, mapply(function(x,y) as.vector(Matrix::crossprod(x, y)), noise_precisions, split(y_minus_field, col(y_minus_field)), SIMPLIFY = F))
  Sys.time()-t1
  
  t1 = Sys.time()
  i_basis_function=1
  eps = parallel::mclapply(
    mc.cores = 5, which(coloring == colo)
    ,
    function(i_basis_function)
    {
      #      # part with observations ####
      B_tau_y_minus_field = Matrix::crossprod(basis_functions[[i_basis_function]], tau_y_minus_field)
      
      B_tau_B_x = list()
      B_tau_B_i = list()
      B_tau_B_j = list()
      for(i_time in seq(time_begin, time_end)){
        B_tau_B = Matrix::crossprod(basis_functions[[i_basis_function]], Matrix::crossprod(noise_precisions[[i_time-time_begin+1]], basis_functions[[i_basis_function]]))
        B_tau_B = Matrix::forceSymmetric(B_tau_B)
        if(length(B_tau_B@x)>0)
        {
          B_tau_B_x = c(B_tau_B_x, list(B_tau_B@x))
          B_tau_B_i = c(B_tau_B_i, list((i_time-time_begin)*ncol(basis_functions[[i_basis_function]]) + B_tau_B@i+1))
          B_tau_B_j = c(B_tau_B_j, list((i_time-time_begin)*ncol(basis_functions[[i_basis_function]]) + findInterval(seq(length(B_tau_B@x))-1,B_tau_B@p[-1])+1))
        }
      }
      B_tau_B_x = unlist(B_tau_B_x)
      B_tau_B_i = unlist(B_tau_B_i)
      B_tau_B_j = unlist(B_tau_B_j)
      
      #Sys.time()-t1
      # part with prior ####
      B_Q_B = lapply(precision_blocks, function(x)Matrix::crossprod(basis_functions[[i_basis_function]], x) %*% basis_functions[[i_basis_function]])
      B_Q_B[[1]]=Matrix::forceSymmetric(B_Q_B[[1]], uplo = "L")
      prior_precision = get_block_toeplitz(B_Q_B, time_end- time_begin+1)
      B_Q_field = (Matrix::crossprod(basis_functions[[i_basis_function]], Q_field[,seq(time_begin, time_end)]))
      
      # posterior ####
      t1 = Sys.time()
      posterior_chol = Matrix::chol(Matrix::sparseMatrix(
        i = c(prior_precision$i, B_tau_B_j), 
        j = c(prior_precision$j, B_tau_B_i), 
        x = c(prior_precision$x, B_tau_B_x), 
        symmetric = T, 
        dims = rep((time_end-time_begin+1)*ncol(basis_functions[[i_basis_function]]), 2)
      ))
      Sys.time()-t1
      
      # sampling coefficients ####
      
      eps = 
        matrix(
          Matrix::solve(
            posterior_chol, 
            rnorm(length(B_tau_y_minus_field))
            + Matrix::solve(
              Matrix::t(posterior_chol), 
              as.vector(
                - B_Q_field 
                - B_tau_y_minus_field
              )
            )
          ), 
          nrow = ncol(basis_functions[[i_basis_function]])
        )
      
      return(eps)
    }, mc.preschedule = F
  )
  Sys.time()-t1
  
  mcmc_nngp_list$chains$chain_1$params$field[,time_begin:time_end] = 
    mcmc_nngp_list$chains$chain_1$params$field[,time_begin:time_end] +
    as.matrix(Reduce("+", mapply(function(e, b) b%*%e, eps, basis_functions[which(coloring==colo)])))
  
  
}



