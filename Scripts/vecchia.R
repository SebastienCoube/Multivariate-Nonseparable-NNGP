


Rcpp::sourceCpp("Scripts/multiply.cpp")
source("Scripts/GMA.R")
source("Scripts/grouping.R")
source("Scripts/multivariate_NN.R")


set.seed(1)

n_loc = 15
n_var = 3
n_time = 3

# spatial locations and variable index
locs = cbind(runif(n_loc), runif(n_loc))
row.names(locs) = seq(n_loc)
var_tag = ceiling(n_var*runif(n_loc));var_tag = match(var_tag, unique(var_tag))

# covariance params
rho = GpGp::exponential_isotropic(c(1, 1, 0), .1*matrix(rnorm(2*n_var), n_var))
rho_vec = rho[lower.tri(rho, diag = F)];  remove(rho)
a2_vec = runif(n_var)
nu_vec = .5 + 2*runif(n_var)
alpha = .0001 * runif(1)
a = runif(1) 
b = runif(1) 
cc = .1 * runif(1) 
lambda = runif(1) 
delta = runif(1) 
r = .1 * runif(1) 
A_vec = runif(n_var)
u = seq(0, n_time-1)

# creating the DAG
DAG = find_ordered_nn_multi(locs = locs, 
                            var_tag = var_tag, 
                            m_whatever_closest = 5, m_same_var = 0, m_other_vars = 0, 
                            lonlat = F)
DAG = guinness_groups(DAG)
summary(sapply(DAG$children, length))
hist(sapply(DAG$children, length), breaks = seq(1, max(sapply(DAG$children, length))))
summary(sapply(DAG$parents, length))


# Pre-computing elements
multiplier = get_multiplier(
  a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
  r = r, A_vec = A_vec, nu_vec = nu_vec, a2_vec = a2_vec, u = u
)
effective_range = get_effective_range(
  a = a, b = b, cc = cc, delta = delta, lambda = lambda, 
  r = r, A_vec = A_vec, a2_vec = a2_vec, u = u
)



my_image = function(m)image(t(m)[,nrow(m):1])

i=1

multi_Vecchia = function(
  DAG, 
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
    #putting parents first, children then so that children arrive at the end of inverse Cholesky factor
    locs_ = locs[c(DAG$parents[[i]],DAG$children[[i]]),] # subset locs
    var_tag_ = var_tag[c(DAG$parents[[i]], DAG$children[[i]])] # subset var tag
    n_loc_ = nrow(locs_)
    
    # computing covmat
    block_lower_tri_idx = lower_tri_idx(n_loc_, diag = T)
    covmat_coeffs = GMA_compressed(
      locs = locs_, 
      var_tag = var_tag_,
      multiplier = multiplier, 
      effective_range = effective_range,
      nu_ = nu_, 
      rho_vec_with_ones = rho_vec_with_ones, 
      n_var = n_var
    )$covmat
    
    covmat_chol = chol(
      expand_block_toeplitz_covmat(
        covmat_coeffs = covmat_coeffs, 
        block_lower_tri_idx = block_lower_tri_idx, 
        n_loc = n_loc_
      )
    )
    
    M = matrix(0, n_time* n_loc_, length(DAG$children[[i]]))
    M[
      cbind(
        seq(nrow(M)-length(DAG$children[[i]])+1, nrow(M)), 
        seq(length(DAG$children[[i]])))
    ]=1
    res$coeffs[[i]] = backsolve(
      r =  covmat_chol, 
      M, transpose = F
    )
  }
  #Sys.time()-t1
  res
}


Linv_coeffs = multi_Vecchia(DAG = DAG, locs = locs, var_tag = var_tag, multiplier = multiplier, effective_range = effective_range, rho_vec = rho_vec, compute_derivative_wrt_rho = F)

put_vecchia_chol_bricks_together = function(Linv_coeffs, DAG, n_time)
{
  mat_idx = list()
  mat_coeffs = list(); for(i in seq(n_time))mat_coeffs[[i]]=list(); names(mat_coeffs)= paste("lag", seq(n_time)-1, sep = "")
  for(i in seq(length(DAG[[1]])))
  {
    n_loc_ = length(DAG$children[[i]])+length(DAG$parents[[i]])
    mat_idx[[i]] = expand.grid(c(DAG$parents[[i]], DAG$children[[i]]), DAG$children[[i]])
    for(j in seq(n_time))
    {
      mat_coeffs[[j]][[i]]= Linv_coeffs$coeffs[[i]][(n_loc_*(n_time-j+1)):(n_loc_*(n_time-j)+1),]  
    }
  }
  row_idx = unlist(lapply(mat_idx, function(x)x[,1]))
  col_idx = unlist(lapply(mat_idx, function(x)x[,2]))
  res = lapply(mat_coeffs, function(x)Matrix::sparseMatrix(i = row_idx, j = col_idx, x = unlist(x)))
  my_image(as.matrix(res[[1]])!=0)
  my_image(as.matrix(res[[2]])!=0)
  my_image(as.matrix(Linv_coeffs$coeffs[[1]])!=0)
}



tatata = Matrix::bdiag(lapply(seq(10000), function(i)matrix(1)))

tatato = Matrix::sparseMatrix(i = ceiling(10000*runif(100000)), 
                              j = ceiling(10000*runif(100000)), 
                              x= rnorm(100000)
                              )
tatato = Matrix::crossprod(tatato)

  tatato = Matrix::sparseMatrix(i = ceiling(10000*runif(100000)), 
                                j = ceiling(10000*runif(100000)), 
                                x= rnorm(100000)
  )
  tatato = Matrix::crossprod(tatato)
t1 = Sys.time()
sapply(seq(1000), function(i)
{
#  tatato+tatata
  tatato[seq(9000, 10000), seq(9000, 10000)]
return(0)
}
  )
Sys.time()-t1

