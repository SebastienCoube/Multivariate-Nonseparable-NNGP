##Rcpp::sourceCpp("Scripts/multiply.cpp")
##source("Scripts/GMA.R")
##source("Scripts/grouping.R")
##source("Scripts/multivariate_NN.R")


# flip image function to match matrix
my_image = function(m)image(t(m)[,nrow(m):1])


# simulate data
set.seed(2)

n_var = 1
n_time = 5

# spatial locations and variable index
locs = as.matrix(expand.grid(seq(40), seq(40)))/20#cbind(runif(n_loc), runif(n_loc))
locs = locs[sample(seq(nrow(locs)), size = nrow(locs), F),]
n_loc =nrow(locs)
row.names(locs) = seq(n_loc)
var_tag = ceiling(n_var*runif(n_loc));var_tag = match(var_tag, unique(var_tag))

# covariance params
rho = GpGp::exponential_isotropic(c(1, 1, 0), .1*matrix(rnorm(2*n_var), n_var))
rho[]=1
rho_vec = rho[lower.tri(rho, diag = F)];  remove(rho)
a2_vec = rep(1, n_var)#.000001*runif(n_var)
nu_vec = rep(1, n_var)#.5 + 2*runif(n_var)
alpha = .01 * runif(1)
a = runif(1) 
b = runif(1) 
cc = 1*runif(1) 
lambda = runif(1) 
delta = runif(1) 
r = 1 * runif(1) 
A_vec = rep(0, n_var)#runif(n_var)
u = seq(0, n_time-1)

# creating the DAG
NNarray_same_time = find_ordered_nn_multi(locs = locs, 
                            var_tag = var_tag, 
                            m_whatever_closest = 10, m_same_var = 5, m_other_vars = 0, 
                            lonlat = F)
NNarray_pevious_times = find_unordered_nn_multi(locs = locs, 
                            var_tag = var_tag, 
                            m_whatever_closest = 10, m_same_var = 5, m_other_vars = 0, 
                            lonlat = F)
DAG = list(
  "children" = lapply(seq(n_loc), function(x)x), 
  "parents_same_time" = NNarray_same_time, 
  "parents_previous_times" = NNarray_pevious_times 
           )

###DAG$children[[11]] = unique(c(DAG$children[[11]], DAG$children[[12]], DAG$children[[13]], DAG$children[[14]]))
###DAG$children[[12]]=NULL; DAG$children[[13]] = NULL; DAG$children[[14]] = NULL
###DAG$parents_same_time[[11]] = unique(c(DAG$parents_same_time[[11]], DAG$parents_same_time[[12]], DAG$parents_same_time[[13]], DAG$parents_same_time[[14]]))
###DAG$parents_same_time[[12]]=NULL; DAG$parents_same_time[[13]] = NULL; DAG$parents_same_time[[14]] = NULL
###DAG$parents_same_time[[11]]=setdiff(DAG$parents_same_time[[11]], DAG$children[[11]])
###DAG$parents_previous_times[[11]] = unique(c(DAG$parents_previous_times[[11]], DAG$parents_previous_times[[12]], DAG$parents_previous_times[[13]], DAG$parents_previous_times[[14]]))
###DAG$parents_previous_times[[12]]=NULL; DAG$parents_previous_times[[13]] = NULL; DAG$parents_previous_times[[14]] = NULL


# getting parts of lower-triangular covmat affected by rho #############
get_rho_idx_from_dag = function(DAG, var_tag)
{
  n_var = max(var_tag)
  idx_1 = mapply(c, DAG$children, DAG$parents_same_time)
  var_idx_current_time = lapply(idx_1, function(idx)position_in_lower_tri_cross_vec(var_tag[idx], n_var = n_var, diag = T))
  var_idx_previous_times = lapply(DAG$parents_previous_times, function(idx)position_in_lower_tri_cross_vec(var_tag[idx], n_var = n_var, diag = F))
  var_idx_cross = mapply(
    function(x, y)outer(x, y, position_in_lower_tri, n_loc = n_var), 
    var_idx_current_time,
    var_idx_previous_times
  )
  
  rho_indices = position_in_lower_tri_cross_vec(seq(n_var), n_var = n_var, diag = F)
  
  res = lapply(rho_indices, function(i)
    list(
      current_time   = lapply(var_idx_current_time,   function(x)get_lower_tri_idx(length(x), T)[which(x==i),]),
      previous_times = lapply(var_idx_previous_times, function(x)get_lower_tri_idx(length(x), T)[which(x==i),]), 
      cross = lapply(var_idx_cross, function(x)cbind(row(x)[x==i], col(x)[x==i]))
    )
    )
  names(res)= outer(seq(n_var), seq(n_var), paste, sep = "_")[lower.tri(matrix(0, n_var, n_var))]
  res
}

rho_idx = get_rho_idx_from_dag(DAG = DAG, var_tag = var_tag )

# get lower triangular indices in the dag

get_lower_tri_idx_DAG = function(DAG)
{
  lower_tri_idx_DAG = list()
  lower_tri_idx_DAG$same_time = lapply(mapply(c,DAG$children, DAG$parents_same_time, SIMPLIFY = F), function(x)get_lower_tri_idx(length(x), T))
  lower_tri_idx_DAG$previous_times = lapply(DAG$parents_previous_times, function(x)get_lower_tri_idx(length(x), T))
  lower_tri_idx_DAG
}

lower_tri_idx_DAG = get_lower_tri_idx_DAG(DAG)

# get var idx in lower tri mat

get_var_idx = function(DAG, var_tag)
{
  n_var = max(var_tag)
  var_idx = list()
  var_idx$current_time = lapply(mapply(c,DAG$children, DAG$parents_same_time, SIMPLIFY = F), function(idx)position_in_lower_tri_cross_vec(var_tag[idx], n_var = n_var, diag = T))
  var_idx$previous_times = lapply(DAG$parents_previous_times, function(idx)position_in_lower_tri_cross_vec(var_tag[idx], n_var = n_var, diag = T))
  var_idx$cross = mapply(
    function(x, y)
    {
      tatato = expand.grid(x, y)
      var_idx_12 = position_in_lower_tri(tatato[,1], tatato[,2], n_var)
    }
      ,
    lapply(mapply(c,DAG$children, DAG$parents_same_time, SIMPLIFY = F), function(x)var_tag[x]),
    lapply(DAG$parents_previous_times                                 , function(x)var_tag[x])
  )
  var_idx
}

var_idx = get_var_idx(DAG, var_tag = var_tag)


# Pre-computing range and multiplier
multiplier = get_multiplier(
  a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
  r = r, A_vec = A_vec, nu_vec = nu_vec, a2_vec = a2_vec, u = u
)
effective_range = get_effective_range(
  a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
  r = r, A_vec = A_vec, a2_vec = a2_vec, u = u
)

##################

i=10
get_linv_coeffs = function(
    DAG, 
    rho_idx, lower_tri_idx_DAG, var_idx,
    locs, var_tag, 
    multiplier, effective_range, 
    rho_vec, 
    compute_derivative_wrt_rho = F
)
{
  res = list()
  res$coeffs = list()
  
  n_time = nrow(multiplier)
  n_var = max(var_tag)
  nu_ = expand_nu(nu_vec)
  rho_vec_with_ones = put_ones_in_rho_vec(rho_vec)
  # looping over each block of the DAG
  #t1 = Sys.time()
  for(i in seq(length(DAG[[1]])))
  {
    print(i)
    covmat_coeffs_same_time = GMA_compressed(
      locs = locs[c(DAG$children[[i]], DAG$parents_same_time[[i]]),, drop = F],
      lower_tri_idx = lower_tri_idx_DAG$same_time[[i]],
      var_idx = var_idx$current_time[[i]],
      multiplier = multiplier[1,,drop=FALSE], 
      effective_range = effective_range,
      nu_ = nu_,
      rho_vec_with_ones = rho_vec_with_ones
    )
    covmat_coeffs_previous = GMA_compressed(
      locs = locs[DAG$parents_previous_times[[i]],, drop = F],
      lower_tri_idx = lower_tri_idx_DAG$previous_times[[i]],
      var_idx = var_idx$previous_times[[i]],
      multiplier = multiplier[-nrow(multiplier),,drop=FALSE], 
      effective_range = effective_range,
      nu_ = nu_,
      rho_vec_with_ones = rho_vec_with_ones
    )
    covmat_coeffs_cross =  
      GMA_rectangular(
        locs_1 = locs[c(DAG$children[[i]], DAG$parents_same_time[[i]]),, drop = F], 
        locs_2 = locs[DAG$parents_previous_times[[i]],, drop = F], 
        var_idx = var_idx$cross[[i]],
        multiplier = multiplier[-1,,drop=FALSE], 
        effective_range = effective_range[-1,,drop=FALSE],
        nu_ = expand_nu(nu_vec),
        n_var =  n_var, 
        rho_vec_with_ones = rho_vec_with_ones
      )
  #my_image( expand_full_covmat(
  #  covmat_previous_periods = expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_previous$covmat,  block_lower_tri_idx = lower_tri_idx_DAG$previous_times[[i]], n_loc = length(DAG$parents_previous_times[[i]])),
  #  covmat_current_period   = expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_same_time$covmat, block_lower_tri_idx = lower_tri_idx_DAG$same_time     [[i]], n_loc = length(DAG$children[[i]])+length(DAG$parents_same_time[[i]])),
  #  side_blocks_rectangles  = covmat_coeffs_cross$covmat
  #))
    covmat_chol = chol(
      expand_full_covmat(
        covmat_previous_periods = diag(nrow(expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_previous$covmat,  block_lower_tri_idx = lower_tri_idx_DAG$previous_times[[i]], n_loc = length(DAG$parents_previous_times[[i]])))),
        covmat_current_period   = expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_same_time$covmat, block_lower_tri_idx = lower_tri_idx_DAG$same_time     [[i]], n_loc = length(DAG$children[[i]])+length(DAG$parents_same_time[[i]])),
        side_blocks_rectangles  = 0*covmat_coeffs_cross$covmat
      )
    )
    M = matrix(0, nrow(covmat_chol), length(DAG$children[[i]]))
    M[cbind(seq(nrow(covmat_chol)- length(DAG$children[[i]])+1, nrow(covmat_chol)), seq(length(DAG$children[[i]])))]=1
    #my_image(covmat_chol)
    res$coeffs[[i]] = backsolve(
      r =  covmat_chol, 
      M
    )
    #plot(res$coeffs[[i]], solve(covmat_chol)[1,])
    # each row of solve covmat_chol factorizes the joint distribuion of children and parents
    # top rows = involves children + parents 
  }
  
  #Sys.time()-t1
  res
}


Linv_coeffs = get_linv_coeffs(
  DAG = DAG, locs = locs, 
  rho_idx = rho_idx, 
  lower_tri_idx_DAG = lower_tri_idx_DAG, 
  var_idx = var_idx, 
  var_tag = var_tag, 
  multiplier = multiplier, 
  effective_range = effective_range, 
  rho_vec = rho_vec, 
  compute_derivative_wrt_rho = F)

coeffs = Linv_coeffs$coeffs

get_vecchia_blocks = function(DAG, coeffs, n_time)
{
  same_time_indices = mapply(c, DAG$parents_same_time, DAG$children, SIMPLIFY = F)
  j = unlist(mapply(function(x, y)outer(x, rep(1, length(y))), same_time_indices, DAG$children))
  i = unlist(mapply(function(x, y)outer(rep(1, length(x)), y), same_time_indices, DAG$children))
  x = unlist(mapply(function(x, y)x[seq(nrow(x)-length(y)+1, nrow(x)),], coeffs, same_time_indices))# selecting bottom of  coeff matrices
  triangular_on_diag = Matrix::sparseMatrix(
    i = i[i>=j],
    j = j[i>=j],
    x = x[i>=j],
    triangular = T
  )
  #j_previous = lapply(
  #  DAG$parents_previous_times, 
  #  function(x)
  #  {
  #    rep(x, n_time-1) + # position within time
  #      max(unlist(DAG$children)) * rep(seq(0, n_time-2), each = length(x)) # switching with respect to time
  #  }
  #)
  #rectangular_below_diag = Matrix::sparseMatrix(
  #  j = unlist(mapply(function(x, y)outer(x, rep(1, length(y))), j_previous, DAG$children)), 
  #  i = unlist(mapply(function(x, y)outer(rep(1, length(x)), y), j_previous, DAG$children)),
  #  x = unlist(mapply(function(x, y)x[seq(nrow(x)-length(y)),], coeffs, same_time_indices, SIMPLIFY = F)) 
  #)
  rectangular_below_diag = 
    lapply(seq(n_time-1), function(idx)
        Matrix::sparseMatrix(
        j = unlist(mapply(function(x, y)outer(x, rep(1, length(y))), DAG$parents_previous_times, DAG$children)), 
        i = unlist(mapply(function(x, y)outer(rep(1, length(x)), y), DAG$parents_previous_times, DAG$children)),
        x = unlist(mapply(function(x, y)x[((idx-1)*length(y)+1):((idx)*length(y)),], coeffs, DAG$parents_previous_times, SIMPLIFY = F), recursive = T) 
      )
      )
  list(
    triangular_on_diag = triangular_on_diag,
    rectangular_below_diag = rectangular_below_diag
  )
  
}  

vecchia_blocks = get_vecchia_blocks(DAG, coeffs, n_time)  


#vecchia_blocks$triangular_on_diag = t(solve(chol(GpGp::exponential_isotropic(c(1,1,0), locs))))
#vecchia_blocks$rectangular_below_diag = lapply(vecchia_blocks$rectangular_below_diag, function(x)0*x)

#my_image(as.matrix(vecchia_blocks$triangular_on_diag)!=0)
#my_image(as.matrix(vecchia_blocks$rectangular_below_diag[[1]])!=0)
#my_image(as.matrix(vecchia_blocks$rectangular_below_diag[[2]])!=0)
  

x = rnorm(100*n_loc)
# multiplicatiion of a vector by the compressed repetitive sparse prior Chol
vecchia_blocks_mult = function(x, vecchia_blocks, n_time)
{
  # putting vector in matrix
  x_ = matrix(x, nrow(vecchia_blocks$triangular_on_diag))
  # result matrix
  res = matrix(0, nrow(x_), ncol(x_))
  # first columns
  res[,1:(n_time-1)]=x_[,1:(n_time-1)]
  # multiplying current time blocks
  res[,n_time:ncol(res)] = as.matrix(vecchia_blocks$triangular_on_diag %*% x_[,n_time:ncol(res)])
  # multiplying previous time blocks
  for(i in seq(n_time-1)) res[,n_time:ncol(res)] = res[,n_time:ncol(res)] + as.matrix(vecchia_blocks$rectangular_below_diag[[i]] %*% x_[,seq(0, ncol(res)-n_time)+i])
  c(res)
}
vecchia_blocks_mult(x, vecchia_blocks, n_time)[-seq(n_time*n_loc)]  -
as.vector(Matrix::bdiag(lapply(seq(length(x)/n_loc), function(i)vecchia_blocks$triangular_on_diag))%*%x)[-seq(n_time*n_loc)]

# solving of a linear by the compressed repetitive sparse prior Chol
vecchia_blocks_solve = function(x, vecchia_blocks, n_time)
{
  # stacking x in  matrix
  x_ = matrix(x, nrow(vecchia_blocks$triangular_on_diag))
  # result
  res = matrix(0, nrow(x_), ncol(x_))
  # initializing result
  res[,1:(n_time-1)]=x_[,1:(n_time-1)]
  # solving recursively
  for(i in seq(n_time, ncol(x_)))
  {
    res[,i] = x_[,i]
    # triangular block matrix solve
    for(j in seq(n_time-1)) res[,i] = res[,i] - 
        as.matrix(vecchia_blocks$rectangular_below_diag[[j]] %*% res[,i-n_time+j])
    res[,i] = as.vector(Matrix::solve(vecchia_blocks$triangular_on_diag, res[,i]))
  }
  c(res)
}

#plot(x, 
#vecchia_blocks_solve(vecchia_blocks_mult(x, vecchia_blocks, n_time), vecchia_blocks, n_time)
#     )
#
#plot(vecchia_blocks_solve(x, vecchia_blocks, n_time))

x = rnorm(100*n_loc)

Bidart::plot_pointillist_painting(locs, vecchia_blocks_solve(x, vecchia_blocks, n_time)[seq(length(x)-n_loc+1, length(x))])



#Bidart::plot_pointillist_painting(locs, as.vector(Matrix::solve(block_diag, rnorm(1000))))
#
#head(block_diag)



