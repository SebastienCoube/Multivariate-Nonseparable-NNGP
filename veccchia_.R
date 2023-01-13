#Rcpp::sourceCpp("Scripts/multiply.cpp")
#source("Scripts/GMA.R")
#source("Scripts/grouping.R")
#source("Scripts/multivariate_NN.R")

# flip image function to match matrix
my_image = function(m)image(t(m)[,nrow(m):1])


# simulate data
set.seed(2)

n_var = 8
n_time = 4

# spatial locations and variable index
locs = as.matrix(expand.grid(seq(100), seq(100)))/20#cbind(runif(n_loc), runif(n_loc))
locs = locs + rnorm(length(locs), 0, .0001)
locs = locs[sample(seq(nrow(locs)), size = nrow(locs), F),]
n_loc =nrow(locs)
row.names(locs) = seq(n_loc)
var_tag = ceiling(n_var*runif(n_loc));var_tag = match(var_tag, unique(var_tag))

# covariance params
rho = GpGp::exponential_isotropic(c(1, 1, 0), matrix(rnorm(2*n_var), n_var))
rho_vec = rho[lower.tri(rho, diag = F)];  remove(rho)
a2_vec = rep(100, n_var)#.000001*runif(n_var)
nu_vec = .5 + 2*runif(n_var)
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


#DAG$children[[10]] = c(10, 11, 12, 13)

# get lower triangular indices in the dag

get_lower_tri_idx_DAG = function(DAG)
{
  lower_tri_idx_DAG = list()
  lower_tri_idx_DAG$same_time = lapply(mapply(c, DAG$parents_same_time, DAG$children, SIMPLIFY = F), function(x)get_lower_tri_idx(length(x), T))
  lower_tri_idx_DAG$previous_times = lapply(DAG$parents_previous_times, function(x)get_lower_tri_idx(length(x), T))
  lower_tri_idx_DAG
}

lower_tri_idx_DAG = get_lower_tri_idx_DAG(DAG)

# get var idx in lower tri mat
get_var_idx = function(DAG, var_tag)
{
  n_var = max(var_tag)
  var_idx = list()
  var_idx$current_time = lapply(mapply(c,DAG$parents_same_time, DAG$children, SIMPLIFY = F), function(idx)position_in_lower_tri_cross_vec(var_tag[idx], n_var = n_var, diag = T))
  var_idx$previous_times = lapply(DAG$parents_previous_times, function(idx)position_in_lower_tri_cross_vec(var_tag[idx], n_var = n_var, diag = T))
  var_idx$cross = mapply(
    function(x, y)
    {
      tatato = expand.grid(x, y)
      var_idx_12 = position_in_lower_tri(tatato[,1], tatato[,2], n_var)
    }
      ,
    lapply(mapply(c, DAG$parents_same_time, DAG$children, SIMPLIFY = F), function(x)var_tag[x]), 
    lapply(DAG$parents_previous_times                                 , function(x)var_tag[x]),# previous time first
    SIMPLIFY = F
  )
  var_idx
}

var_idx = get_var_idx(DAG, var_tag = var_tag)


# get indices in covariance matrix that are touched by a move in rho

get_rho_idx = function(DAG, var_idx, lower_tri_idx_DAG, n_var)
{
  # adding rectangular indices to lower tri indices
  indices = 
    c(lower_tri_idx_DAG, 
      "cross" = list(mapply(
        function(x, y)as.matrix(expand.grid(seq(length(x)), seq(length(y)))),
        mapply(c, DAG$parents_same_time, DAG$children), 
        DAG$parents_previous_times,
        SIMPLIFY = F
      ))
    )
  # get combination idx of the variable pairs in the indices arrays
##  combination_idx = 
##    lapply(indices, function(idx_list)
##      lapply(idx_list, function(idx_array)
##        position_in_lower_tri(var_tag[idx_array[,1]], var_tag[idx_array[,2]], n_var)
##      )
##    )
  res = lapply(seq(n_var*(n_var+1)/2)[-cumsum(c(1, seq(n_var, 2)))], # indices of rho after diagonal terms are removed
               function(i)
                 mapply(indices, var_idx, FUN = function(list1, list2)
                   mapply(list1, list2, FUN = function(idx_array, var_combination_idx)idx_array[var_combination_idx==i,,drop=F], SIMPLIFY = F), 
                   SIMPLIFY = F
                   )
  )
  names(res)=paste("rho", seq(n_var*(n_var+1)/2)[-cumsum(c(1, seq(n_var, 2)))], sep = "_")
  res
}

rho_idx = get_rho_idx(DAG = DAG, var_idx = var_idx, lower_tri_idx_DAG = lower_tri_idx_DAG, n_var = n_var)

#### # Check that the right indices are selected by rho_idx
####plot(indices$cross[[100]], col = var_idx$cross[[100]])
####points(rho_idx$rho_2$cross[[100]], pch = 3)
####
####plot(indices$same_time[[100]], col = var_idx$current_time[[100]])
####points(rho_idx$rho_2$same_time[[100]], pch = 3)

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
    locs, var_tag, 
    DAG, 
    lower_tri_idx_DAG, 
    var_idx,
    multiplier, effective_range, 
    rho_vec
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
  res = parallel::mclapply(
    X = seq(length(DAG[[1]])), mc.cores = parallel::detectCores(), 
    FUN = function(i)
    {
      #print(i)
      covmat_coeffs_same_time = GMA_compressed(
        locs = locs[c(DAG$parents_same_time[[i]], DAG$children[[i]]),, drop = F],
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
          locs_1 = locs[c(DAG$parents_same_time[[i]], DAG$children[[i]]),, drop = F], 
          locs_2 = locs[DAG$parents_previous_times[[i]],, drop = F], 
          var_idx = var_idx$cross[[i]],
          multiplier = multiplier[-1,,drop=FALSE], 
          effective_range = effective_range[-1,,drop=FALSE],
          nu_ = expand_nu(nu_vec),
          n_var =  n_var, 
          rho_vec_with_ones = rho_vec_with_ones
        )
      # my_image( expand_full_covmat(
      #   covmat_previous_periods = expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_previous$covmat,  block_lower_tri_idx = lower_tri_idx_DAG$previous_times[[i]], n_loc = length(DAG$parents_previous_times[[i]])),
      #   covmat_current_period   = expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_same_time$covmat, block_lower_tri_idx = lower_tri_idx_DAG$same_time     [[i]], n_loc = length(DAG$children[[i]])+length(DAG$parents_same_time[[i]])),
      #   side_blocks_rectangles  = covmat_coeffs_cross$covmat
      # ))
      covmat_chol = chol(
        expand_full_covmat(
          covmat_previous_periods = expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_previous$covmat,  block_lower_tri_idx = lower_tri_idx_DAG$previous_times[[i]], n_loc = length(DAG$parents_previous_times[[i]])),
          covmat_current_period   = expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_same_time$covmat, block_lower_tri_idx = lower_tri_idx_DAG$same_time     [[i]], n_loc = length(DAG$children[[i]])+length(DAG$parents_same_time[[i]])),
          side_blocks_rectangles  = covmat_coeffs_cross$covmat
        )
      )
      M = matrix(0, nrow(covmat_chol), length(DAG$children[[i]]))
      M[cbind(seq(nrow(covmat_chol)- length(DAG$children[[i]])+1, nrow(covmat_chol)), seq(length(DAG$children[[i]])))]=1
      return(backsolve(
        r =  covmat_chol, 
        M
      ))
      #my_image(covmat_chol)
      ### full_covmat = expand_full_covmat(
      ###   covmat_previous_periods = expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_previous$covmat,  block_lower_tri_idx = lower_tri_idx_DAG$previous_times[[i]], n_loc = length(DAG$parents_previous_times[[i]])),
      ###   covmat_current_period   = expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_same_time$covmat, block_lower_tri_idx = lower_tri_idx_DAG$same_time     [[i]], n_loc = length(DAG$children[[i]])+length(DAG$parents_same_time[[i]])),
      ###   side_blocks_rectangles  = 0*covmat_coeffs_cross$covmat
      ### )
      ### 1/sqrt(full_covmat[nrow(full_covmat), nrow(full_covmat)] - full_covmat[nrow(full_covmat), -nrow(full_covmat)]%*%solve(full_covmat[-nrow(full_covmat), -nrow(full_covmat)]) %*%full_covmat[nrow(full_covmat), -nrow(full_covmat)])
      ### (solve(full_covmat[-nrow(full_covmat), -nrow(full_covmat)]) %*%full_covmat[nrow(full_covmat), -nrow(full_covmat)])/c(sqrt(full_covmat[nrow(full_covmat), nrow(full_covmat)] - full_covmat[nrow(full_covmat), -nrow(full_covmat)]%*%solve(full_covmat[-nrow(full_covmat), -nrow(full_covmat)]) %*%full_covmat[nrow(full_covmat), -nrow(full_covmat)]))
      
      #plot(res$coeffs[[i]], solve(covmat_chol)[1,])
      # each row of solve covmat_chol factorizes the joint distribuion of children and parents
      # top rows = involves children + parents 
    }
    )
  #Sys.time()-t1
  res
}


get_sigma_12 = function(
      coeffs_same_time, 
      block_lower_tri_idx_same_time, n_loc_same_time, 
      coeffs_cross, n_children
){
  res = matrix(0, dim(coeffs_cross)[2]*dim(coeffs_cross)[3]+n_loc_same_time, n_children)
  res[seq(dim(coeffs_cross)[3]*dim(coeffs_cross)[2]), ] = coeffs_cross[(dim(coeffs_cross)[1]-n_children+1):dim(coeffs_cross)[1],,]
  res[-seq(dim(coeffs_cross)[3]*dim(coeffs_cross)[2]), ] = 
    expand_covmat_into_blocks(
      coeffs_same_time, 
      block_lower_tri_idx = block_lower_tri_idx_same_time,
      n_loc = n_loc_same_time)[[1]][,(n_loc_same_time-n_children+1):n_loc_same_time]
  res
}

x_plus_tx = function(x)x+t(x)

i=10
get_linv_coeffs_and_derivative = function(
    locs, var_tag, 
    DAG, 
    rho_idx, lower_tri_idx_DAG, var_idx,
    multiplier, effective_range, 
    rho_vec, 
    cl)
{
  n_time = nrow(multiplier)
  n_var = max(var_tag)
  nu_ = expand_nu(nu_vec)
  rho_vec_with_ones = put_ones_in_rho_vec(rho_vec)
  # looping over each block of the DAG
  #t1 = Sys.time()
  
  final_res = parallel::mclapply(seq(length(DAG[[1]])), mc.cores = parallel::detectCores(),
  FUN = function(i)
  {
    res = list()
    n_children = length(DAG$children[[i]])
    n_loc_same_time = length(DAG$children[[i]])+length(DAG$parents_same_time[[i]])
    n_loc_previous = length(DAG$parents_previous_times[[i]])
    # computing covmat coefficients
    covmat_coeffs_same_time = GMA_compressed(
      locs = locs[c(DAG$parents_same_time[[i]], DAG$children[[i]]),, drop = F],
      lower_tri_idx = lower_tri_idx_DAG$same_time[[i]],
      var_idx = var_idx$current_time[[i]],
      multiplier = multiplier[1,,drop=FALSE], 
      effective_range = effective_range,
      nu_ = nu_,
      rho_vec_with_ones = rho_vec_with_ones, return_with_rho = T
    )
    covmat_coeffs_previous = GMA_compressed(
      locs = locs[DAG$parents_previous_times[[i]],, drop = F],
      lower_tri_idx = lower_tri_idx_DAG$previous_times[[i]],
      var_idx = var_idx$previous_times[[i]],
      multiplier = multiplier[-nrow(multiplier),,drop=FALSE], 
      effective_range = effective_range,
      nu_ = nu_,
      rho_vec_with_ones = rho_vec_with_ones, return_with_rho = T
    )
    covmat_coeffs_cross =  
      GMA_rectangular(
        locs_1 = locs[c(DAG$parents_same_time[[i]], DAG$children[[i]]),, drop = F], 
        locs_2 = locs[DAG$parents_previous_times[[i]],, drop = F], 
        var_idx = var_idx$cross[[i]],
        multiplier = multiplier[-1,,drop=FALSE], 
        effective_range = effective_range[-1,,drop=FALSE],
        nu_ = expand_nu(nu_vec),
        n_var =  n_var, 
        rho_vec_with_ones = rho_vec_with_ones, 
        return_with_rho = T
      )
    # my_image( expand_full_covmat(
    #   covmat_previous_periods = expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_previous$covmat,  block_lower_tri_idx = lower_tri_idx_DAG$previous_times[[i]], n_loc = length(DAG$parents_previous_times[[i]])),
    #   covmat_current_period   = expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_same_time$covmat, block_lower_tri_idx = lower_tri_idx_DAG$same_time     [[i]], n_loc = length(DAG$children[[i]])+length(DAG$parents_same_time[[i]])),
    #   side_blocks_rectangles  = covmat_coeffs_cross$covmat
    # ))
    # computing Cholesky decomposition of covmat
    covmat = expand_full_covmat(
      covmat_previous_periods = expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_previous$covmat,  block_lower_tri_idx = lower_tri_idx_DAG$previous_times[[i]], n_loc = n_loc_previous),
      covmat_current_period   = expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_same_time$covmat, block_lower_tri_idx = lower_tri_idx_DAG$same_time     [[i]], n_loc = n_loc_same_time),
      side_blocks_rectangles  = covmat_coeffs_cross$covmat
    )
    covmat_chol = chol(covmat)
    Sigma_12 = covmat[,seq(nrow(covmat)-n_children+1,nrow(covmat)), drop = F]
    Sigma_12[(nrow(Sigma_12)-n_children+1):nrow(Sigma_12),]=0
    # computing covmat coefficients
    M = matrix(0, nrow(covmat_chol), length(DAG$children[[i]]))
    M[cbind(seq(nrow(covmat_chol)- length(DAG$children[[i]])+1, nrow(covmat_chol)), seq(length(DAG$children[[i]])))]=1
    res$coeffs = backsolve(
      r =  covmat_chol, 
      M
    )
    # precision of the children knowing the parents
    cond_precision_chol = (res$coeffs[-seq(nrow(res$coeffs)-n_children),])
    # solve(Sigma_1_1) Sigma 12, called "salt" because it goes everywhere
    salt = backsolve(
      (covmat_chol), transpose = T,
      Sigma_12
    )
    salt[seq(nrow(salt)-n_children+1, nrow(salt)),]=0
    salt = 
      backsolve(
        (covmat_chol), salt
      )
    salt[seq(nrow(salt)-n_children+1, nrow(salt)),]=0
    # looping on rho
    res$coeffs_derivatives = list()
    for(j in seq(length(rho_idx)))
    {
      # solve(Sigma 11) (d Sigma 11 / d rho) *  solve(Sigma 11) Sigma 12 ########
      spam_sandwich = apply(
        salt, 
        2, 
        function(x)multiply_vector_full_covmat_sparse(
          v = x,
          covmat_coeffs_current_period = covmat_coeffs_same_time$covmat_without_rho, 
          idx_mat_current_period = rho_idx[[j]]$same_time[[i]], 
          n_loc_current_period = n_loc_same_time, 
          covmat_coeffs_previous_period = covmat_coeffs_previous$covmat_without_rho, 
          n_loc_previous_period = n_loc_previous, 
          idx_mat_previous_period = rho_idx[[j]]$previous_times[[i]], 
          side_blocks_rectangles = covmat_coeffs_cross$covmat_without_rho, 
          idx_mat_side_block_rectangles = rho_idx[[j]]$cross[[i]]
        )
      )
      spam_sandwich[seq(nrow(salt)-n_children+1, nrow(salt)),]=0
      spam_sandwich =backsolve(
        (covmat_chol), transpose = T, spam_sandwich
      )
      spam_sandwich[seq(nrow(salt)-n_children+1, nrow(salt)),]=0
      spam_sandwich = 
        backsolve(
        (covmat_chol),spam_sandwich
        )
      spam_sandwich[seq(nrow(salt)-n_children+1, nrow(salt)),]=0
      
      # - solve(Sigma 11) d Sigma 12 / d rho #########
      pickle = 
        # (d Sigma 12 / d rho) is here
        apply(
          M, 
          2, 
          function(x)multiply_vector_full_covmat_sparse(
            v = x,
            covmat_coeffs_current_period = covmat_coeffs_same_time$covmat_without_rho, 
            idx_mat_current_period = rho_idx[[j]]$same_time[[i]], 
            n_loc_current_period = n_loc_same_time, 
            covmat_coeffs_previous_period = NULL, 
            n_loc_previous_period = n_loc_previous, 
            idx_mat_previous_period = NULL, 
            side_blocks_rectangles = covmat_coeffs_cross$covmat_without_rho, 
            idx_mat_side_block_rectangles = rho_idx[[j]]$cross[[i]]
          )
        )
      pickle[seq(nrow(salt)-n_children+1, nrow(salt)),]=0
      pickle=
          backsolve(
            (covmat_chol), transpose = T,
            pickle
          )
      pickle[seq(nrow(salt)-n_children+1, nrow(salt)),]=0
      pickle = 
        backsolve(
          (covmat_chol), pickle
          )
      pickle[seq(nrow(salt)-n_children+1, nrow(salt)),]=0
      ####### derivative of the conitional variance multiplicator ########
      # d(solve(chol(Sigma 22  - Sigma 12 solve(Sigma 11) Sigma 21))) / d rho = 
      # (solve(chol(Sigma 22  - Sigma 12 solve(Sigma 11) Sigma 21 + 
      # eps * d(Sigma 22  - Sigma 12 Sigma 11 Sigma 21) / d rho ))) - 
      # solve(chol(Sigma 22  - Sigma 12 Sigma 11 Sigma 21))) / eps 
      derivative_cond_precision_chol = 10000 *
        (
          t(solve(chol(
            tcrossprod(solve(t(cond_precision_chol)))  # Sigma 22  - Sigma 12 solve(Sigma 11) Sigma 21
            + .0001 * 
              apply( # d Sigma 22  / d rho
                M, 
                2, 
                function(x)multiply_vector_full_covmat_sparse(
                  v = x,
                  covmat_coeffs_current_period = covmat_coeffs_same_time$covmat_without_rho, 
                  idx_mat_current_period = rho_idx[[j]]$same_time[[i]], 
                  n_loc_current_period = n_loc_same_time, 
                  covmat_coeffs_previous_period = NULL, 
                  idx_mat_previous_period = NULL, 
                  idx_mat_side_block_rectangles = NULL
                )
              )[seq(nrow(M)-ncol(M)+1, nrow(M)),]
            - .0001 * x_plus_tx(crossprod(pickle, Sigma_12))# d Sigma 12 / d rho solve(Sigma 11) Sigma 12 + Sigma 12 solve(Sigma 11) d Sigma 12 / d rho
            + .0001 * crossprod(spam_sandwich, Sigma_12) # Sigma 12 d Sigma 11/ d rho  Sigma 12
          )))- 
            t(cond_precision_chol)
        )
      
      ###### putting stuff together ######
      # d(solve(chol(cond var)))
      res$coeffs_derivatives[[j]] = - salt 
      # adding diag at the bottom
      res$coeffs_derivatives[[j]][
        cbind(
        seq(nrow(pickle)-ncol(pickle)+1, nrow(pickle)), 
        seq(ncol(pickle))
        )]  =1
      # multiplying by conditional variance cholesky derivative
      res$coeffs_derivatives[[j]] = res$coeffs_derivatives[[j]] %*% t(derivative_cond_precision_chol)
      # adding derivative of conditional mean multiplicator
      res$coeffs_derivatives[[j]] = res$coeffs_derivatives[[j]] + (spam_sandwich - pickle) %*% (cond_precision_chol)
        ####################
    }
  res
  }
  )
  #Sys.time()-t1
  final_res
}


tatato = get_linv_coeffs_and_derivative(
  DAG = DAG, locs = locs, 
  lower_tri_idx_DAG = lower_tri_idx_DAG, 
  var_idx = var_idx, 
  var_tag = var_tag, 
  multiplier = multiplier, 
  effective_range = effective_range, 
  rho_vec = rho_vec, rho_idx = rho_idx, cl = cl)

tatata = get_linv_coeffs_and_derivative(
  DAG = DAG, locs = locs, 
  lower_tri_idx_DAG = lower_tri_idx_DAG, 
  var_idx = var_idx, 
  var_tag = var_tag, 
  multiplier = multiplier, 
  effective_range = effective_range, 
  rho_vec = rho_vec+.0001, rho_idx = rho_idx)

((tatata[[10]]$coeffs - tatato[[10]]$coeffs)*10000)/
(tatato[[10]]$coeffs_derivatives[[1]])

i=10
get_linv_coeffs_and_derivative_ = function(
    locs, var_tag, 
    DAG, 
    rho_idx, lower_tri_idx_DAG, var_idx,
    multiplier, effective_range, 
    rho_vec, 
    cl)
{
  n_time = nrow(multiplier)
  n_var = max(var_tag)
  nu_ = expand_nu(nu_vec)
  rho_vec_with_ones = put_ones_in_rho_vec(rho_vec)
  # looping over each block of the DAG
  #t1 = Sys.time()
  
  final_res = parallel::mclapply(seq(length(DAG[[1]])), mc.cores = parallel::detectCores(),
  FUN = function(i)
  {
    res = list()
    n_children = length(DAG$children[[i]])
    n_loc_same_time = length(DAG$children[[i]])+length(DAG$parents_same_time[[i]])
    n_loc_previous = length(DAG$parents_previous_times[[i]])
    # computing covmat coefficients
    covmat_coeffs_same_time = GMA_compressed(
      locs = locs[c(DAG$parents_same_time[[i]], DAG$children[[i]]),, drop = F],
      lower_tri_idx = lower_tri_idx_DAG$same_time[[i]],
      var_idx = var_idx$current_time[[i]],
      multiplier = multiplier[1,,drop=FALSE], 
      effective_range = effective_range,
      nu_ = nu_,
      rho_vec_with_ones = rho_vec_with_ones, return_with_rho = T
    )
    covmat_coeffs_previous = GMA_compressed(
      locs = locs[DAG$parents_previous_times[[i]],, drop = F],
      lower_tri_idx = lower_tri_idx_DAG$previous_times[[i]],
      var_idx = var_idx$previous_times[[i]],
      multiplier = multiplier[-nrow(multiplier),,drop=FALSE], 
      effective_range = effective_range,
      nu_ = nu_,
      rho_vec_with_ones = rho_vec_with_ones, return_with_rho = T
    )
    covmat_coeffs_cross =  
      GMA_rectangular(
        locs_1 = locs[c(DAG$parents_same_time[[i]], DAG$children[[i]]),, drop = F], 
        locs_2 = locs[DAG$parents_previous_times[[i]],, drop = F], 
        var_idx = var_idx$cross[[i]],
        multiplier = multiplier[-1,,drop=FALSE], 
        effective_range = effective_range[-1,,drop=FALSE],
        nu_ = expand_nu(nu_vec),
        n_var =  n_var, 
        rho_vec_with_ones = rho_vec_with_ones, 
        return_with_rho = T
      )
    # my_image( expand_full_covmat(
    #   covmat_previous_periods = expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_previous$covmat,  block_lower_tri_idx = lower_tri_idx_DAG$previous_times[[i]], n_loc = length(DAG$parents_previous_times[[i]])),
    #   covmat_current_period   = expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_same_time$covmat, block_lower_tri_idx = lower_tri_idx_DAG$same_time     [[i]], n_loc = length(DAG$children[[i]])+length(DAG$parents_same_time[[i]])),
    #   side_blocks_rectangles  = covmat_coeffs_cross$covmat
    # ))
    # computing Cholesky decomposition of covmat
    covmat = expand_full_covmat(
      covmat_previous_periods = expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_previous$covmat,  block_lower_tri_idx = lower_tri_idx_DAG$previous_times[[i]], n_loc = n_loc_previous),
      covmat_current_period   = expand_block_toeplitz_covmat(covmat_coeffs = covmat_coeffs_same_time$covmat, block_lower_tri_idx = lower_tri_idx_DAG$same_time     [[i]], n_loc = n_loc_same_time),
      side_blocks_rectangles  = covmat_coeffs_cross$covmat
    )
    covmat_chol = chol(covmat)
    Sigma_12 = covmat[,seq(nrow(covmat)-n_children+1,nrow(covmat)), drop = F]
    Sigma_12[(nrow(Sigma_12)-n_children+1):nrow(Sigma_12),]=0
    # computing covmat coefficients
    M = matrix(0, nrow(covmat_chol), length(DAG$children[[i]]))
    M[cbind(seq(nrow(covmat_chol)- length(DAG$children[[i]])+1, nrow(covmat_chol)), seq(length(DAG$children[[i]])))]=1
    res$coeffs = backsolve(
      r =  covmat_chol, 
      M
    )
    # precision of the children knowing the parents
    cond_precision_chol = (res$coeffs[-seq(nrow(res$coeffs)-n_children),])
    # solve(Sigma_1_1) Sigma 12, called "salt" because it goes everywhere
    salt = backsolve(
      (covmat_chol), transpose = T,
      Sigma_12
    )
    salt[seq(nrow(salt)-n_children+1, nrow(salt)),]=0
    salt = 
      backsolve(
        (covmat_chol), salt
      )
    salt[seq(nrow(salt)-n_children+1, nrow(salt)),]=0
    # looping on rho
    res$coeffs_derivatives = list()
##    for(j in seq(length(rho_idx)))
##    {
##      # solve(Sigma 11) (d Sigma 11 / d rho) *  solve(Sigma 11) Sigma 12 ########
##      spam_sandwich = apply(
##        salt, 
##        2, 
##        function(x)multiply_vector_full_covmat_sparse(
##          v = x,
##          covmat_coeffs_current_period = covmat_coeffs_same_time$covmat_without_rho, 
##          idx_mat_current_period = rho_idx[[j]]$same_time[[i]], 
##          n_loc_current_period = n_loc_same_time, 
##          covmat_coeffs_previous_period = covmat_coeffs_previous$covmat_without_rho, 
##          n_loc_previous_period = n_loc_previous, 
##          idx_mat_previous_period = rho_idx[[j]]$previous_times[[i]], 
##          side_blocks_rectangles = covmat_coeffs_cross$covmat_without_rho, 
##          idx_mat_side_block_rectangles = rho_idx[[j]]$cross[[i]]
##        )
##      )
##      spam_sandwich[seq(nrow(salt)-n_children+1, nrow(salt)),]=0
##      spam_sandwich =backsolve(
##        (covmat_chol), transpose = T, spam_sandwich
##      )
##      spam_sandwich[seq(nrow(salt)-n_children+1, nrow(salt)),]=0
##      spam_sandwich = 
##        backsolve(
##        (covmat_chol),spam_sandwich
##        )
##      spam_sandwich[seq(nrow(salt)-n_children+1, nrow(salt)),]=0
##      
##      # - solve(Sigma 11) d Sigma 12 / d rho #########
##      pickle = 
##        # (d Sigma 12 / d rho) is here
##        apply(
##          M, 
##          2, 
##          function(x)multiply_vector_full_covmat_sparse(
##            v = x,
##            covmat_coeffs_current_period = covmat_coeffs_same_time$covmat_without_rho, 
##            idx_mat_current_period = rho_idx[[j]]$same_time[[i]], 
##            n_loc_current_period = n_loc_same_time, 
##            covmat_coeffs_previous_period = NULL, 
##            n_loc_previous_period = n_loc_previous, 
##            idx_mat_previous_period = NULL, 
##            side_blocks_rectangles = covmat_coeffs_cross$covmat_without_rho, 
##            idx_mat_side_block_rectangles = rho_idx[[j]]$cross[[i]]
##          )
##        )
##      pickle[seq(nrow(salt)-n_children+1, nrow(salt)),]=0
##      pickle=
##          backsolve(
##            (covmat_chol), transpose = T,
##            pickle
##          )
##      pickle[seq(nrow(salt)-n_children+1, nrow(salt)),]=0
##      pickle = 
##        backsolve(
##          (covmat_chol), pickle
##          )
##      pickle[seq(nrow(salt)-n_children+1, nrow(salt)),]=0
##      ####### derivative of the conitional variance multiplicator ########
##      # d(solve(chol(Sigma 22  - Sigma 12 solve(Sigma 11) Sigma 21))) / d rho = 
##      # (solve(chol(Sigma 22  - Sigma 12 solve(Sigma 11) Sigma 21 + 
##      # eps * d(Sigma 22  - Sigma 12 Sigma 11 Sigma 21) / d rho ))) - 
##      # solve(chol(Sigma 22  - Sigma 12 Sigma 11 Sigma 21))) / eps 
##      derivative_cond_precision_chol = 10000 *
##        (
##          t(solve(chol(
##            tcrossprod(solve(t(cond_precision_chol)))  # Sigma 22  - Sigma 12 solve(Sigma 11) Sigma 21
##            + .0001 * 
##              apply( # d Sigma 22  / d rho
##                M, 
##                2, 
##                function(x)multiply_vector_full_covmat_sparse(
##                  v = x,
##                  covmat_coeffs_current_period = covmat_coeffs_same_time$covmat_without_rho, 
##                  idx_mat_current_period = rho_idx[[j]]$same_time[[i]], 
##                  n_loc_current_period = n_loc_same_time, 
##                  covmat_coeffs_previous_period = NULL, 
##                  idx_mat_previous_period = NULL, 
##                  idx_mat_side_block_rectangles = NULL
##                )
##              )[seq(nrow(M)-ncol(M)+1, nrow(M)),]
##            - .0001 * x_plus_tx(crossprod(pickle, Sigma_12))# d Sigma 12 / d rho solve(Sigma 11) Sigma 12 + Sigma 12 solve(Sigma 11) d Sigma 12 / d rho
##            + .0001 * crossprod(spam_sandwich, Sigma_12) # Sigma 12 d Sigma 11/ d rho  Sigma 12
##          )))- 
##            t(cond_precision_chol)
##        )
##      
##      ###### putting stuff together ######
##      # d(solve(chol(cond var)))
##      res$coeffs_derivatives[[j]] = - salt 
##      # adding diag at the bottom
##      res$coeffs_derivatives[[j]][
##        cbind(
##        seq(nrow(pickle)-ncol(pickle)+1, nrow(pickle)), 
##        seq(ncol(pickle))
##        )]  =1
##      # multiplying by conditional variance cholesky derivative
##      res$coeffs_derivatives[[j]] = res$coeffs_derivatives[[j]] %*% t(derivative_cond_precision_chol)
##      # adding derivative of conditional mean multiplicator
##      res$coeffs_derivatives[[j]] = res$coeffs_derivatives[[j]] + (spam_sandwich - pickle) %*% (cond_precision_chol)
##        ####################
##    }
  res
  }
  )
  #Sys.time()-t1
  final_res
}


tatato = get_linv_coeffs_and_derivative_(
  DAG = DAG, locs = locs, 
  lower_tri_idx_DAG = lower_tri_idx_DAG, 
  var_idx = var_idx, 
  var_tag = var_tag, 
  multiplier = multiplier, 
  effective_range = effective_range, 
  rho_vec = rho_vec, rho_idx = rho_idx, cl = cl)



##############

Linv_coeffs = get_linv_coeffs(
  DAG = DAG, locs = locs, 
  lower_tri_idx_DAG = lower_tri_idx_DAG, 
  var_idx = var_idx, 
  var_tag = var_tag, 
  multiplier = multiplier, 
  effective_range = effective_range, 
  rho_vec = rho_vec)
 
coeffs = Linv_coeffs

get_vecchia_blocks = function(DAG, coeffs, n_time)
{
  same_time_indices = mapply(c, DAG$parents_same_time, DAG$children, SIMPLIFY = F)
  j = unlist(mapply(function(x, y)outer(x, rep(1, length(y))), same_time_indices, DAG$children))
  i = unlist(mapply(function(x, y)outer(rep(1, length(x)), y), same_time_indices, DAG$children))
  x = unlist(mapply(function(x, y)x[seq(nrow(x)-length(y)+1, nrow(x)),], coeffs, same_time_indices))# selecting bottom of  coeff matrices
  tatato = unlist(mapply(function(x, y)x[-seq(nrow(x)-length(y)+1, nrow(x)),], coeffs, same_time_indices))# selecting bottom of  coeff matrices
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
        x = unlist(mapply(function(x, y)x[seq((idx-1)*length(y)+1,(idx)*length(y)),], coeffs, DAG$parents_previous_times, SIMPLIFY = F), recursive = T) 
      )
      )
  names(rectangular_below_diag)=sapply(seq(n_time-1), function(i)paste("lag=", n_time-i, sep=""))
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

# multiplicatiion of a vector by the compressed repetitive sparse prior Chol
vecchia_blocks_t_mult = function(x, vecchia_blocks, n_time)
{
  # putting vector in matrix
  x_ = matrix(x, nrow(vecchia_blocks$triangular_on_diag))
  # result matrix
  res = matrix(0, nrow(x_), ncol(x_))
  # first columns
  res[,1:(n_time-1)]=x_[,1:(n_time-1)]
  # multiplying current time blocks
  res[,n_time:ncol(res)] = as.matrix(t(x_[,n_time:ncol(res)]) %*% vecchia_blocks$triangular_on_diag)
  # multiplying previous time blocks
  for(i in seq(n_time-1)) res[,n_time:ncol(res)] = res[,n_time:ncol(res)] + t(as.matrix(t(x_[,seq(0, ncol(res)-n_time)+i]) %*% vecchia_blocks$rectangular_below_diag[[i]]))
  c(res)
}


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

# checking t mult and mult
x1 = rnorm(3600*n_loc)
x2 = rnorm(3600*n_loc)
vecchia_blocks_t_mult(x1, vecchia_blocks, n_time)*x2 - 
  vecchia_blocks_mult(x2, vecchia_blocks, n_time)*x1

# checking solve and mult
x = rnorm(3600*n_loc)
x - vecchia_blocks_mult(vecchia_blocks_solve(x, vecchia_blocks, n_time), vecchia_blocks, n_time)





w = rnorm(20*n_loc)
w[seq(10*n_loc)]=0


Linv_coeffs_no_rho = get_linv_coeffs(
  DAG = DAG, locs = locs, 
  lower_tri_idx_DAG = lower_tri_idx_DAG, 
  var_idx = var_idx, 
  var_tag = var_tag, 
  multiplier = multiplier, 
  effective_range = effective_range, 
  rho_vec = rho_vec/rho_vec)

vecchia_blocks_no_rho = get_vecchia_blocks(DAG, Linv_coeffs_no_rho, n_time)

#Need to do solve from bottom with vecchia blocks no rho

derivative_sufficient_dens_wrt_rho = function(
  w, vecchia_blocks, vecchia_blocks_no_rho, DAG, n_time, n_var, var_tag
  )
{
  # d/d rho wT RTR w = d/d rho wT ((RTR)-1)-1 w 
  # = - wT ((RTR)-1)-1 d/d rho (RTR)-1 ((RTR)-1)-1 w
  # ~= - wT ((RTR)-1)-1 d/d rho  (rho_mat o (R0TR0)-1) ((RTR)-1)-1 w
  # = - wT ((RTR)-1)-1   (idx (vat i x var j) o (R0TR0)-1) ((RTR)-1)-1 w
  # = - wT RTR   (idx (vat i x var j) o (R0TR0)-1) RTR w
  # = - (RTR w o 1 var i) ((R0TR0)-1) (RTR w o 1 var j) - (RTR w o 1 var j) ((R0TR0)-1) (RTR w o 1 var i)
  RT_R_w = vecchia_blocks_t_mult(vecchia_blocks_mult(x = w, vecchia_blocks = vecchia_blocks, n_time = n_time), vecchia_blocks = vecchia_blocks, n_time = n_time)
  var_tag_mat = matrix(0, length(RT_R_w), ncol = n_var)
  var_tag_mat[cbind(seq(nrow(var_tag_mat)), rep(var_tag, nrow(var_tag_mat)/length(var_tag)))]=1
  derivative_mat = -.5 * crossprod(apply(var_tag_mat*RT_R_w, 2, function(x)vecchia_blocks_solve(x, vecchia_blocks = vecchia_blocks_no_rho, n_time = n_time)))
}

covmat_no_rho = tcrossprod(apply(diag(20*n_loc), 2, function(x)vecchia_blocks_solve(x, vecchia_blocks = vecchia_blocks_no_rho, n_time = n_time)))
my_image(covmat_no_rho)

(RT_R_w*var_tag_mat[,1]) %*%
covmat_no_rho %*%
(RT_R_w*var_tag_mat[,2])

crossprod(
  vecchia_blocks_solve(RT_R_w*var_tag_mat[,1], vecchia_blocks = vecchia_blocks_no_rho, n_time = n_time),
  vecchia_blocks_solve(RT_R_w*var_tag_mat[,2], vecchia_blocks = vecchia_blocks_no_rho, n_time = n_time)
  )

Linv_coeffs = get_linv_coeffs(
  DAG = DAG, locs = locs, 
  lower_tri_idx_DAG = lower_tri_idx_DAG, 
  var_idx = var_idx, 
  var_tag = var_tag, 
  multiplier = multiplier, 
  effective_range = effective_range, 
  rho_vec = rho_vec)
vecchia_blocks = get_vecchia_blocks(DAG, Linv_coeffs, n_time)  
tatato = rep(0, length(rho_vec))
for(i in seq(length(rho_vec)))
{
  rho_vec_ = rho_vec; rho_vec_[i]=rho_vec_[i]+.001
  Linv_coeffs_ = get_linv_coeffs(
    DAG = DAG, locs = locs, 
    lower_tri_idx_DAG = lower_tri_idx_DAG, 
    var_idx = var_idx, 
    var_tag = var_tag, 
    multiplier = multiplier, 
    effective_range = effective_range, 
    rho_vec = rho_vec_)
  vecchia_blocks_ = get_vecchia_blocks(DAG, Linv_coeffs_, n_time)  
  tatato[i] = 
  1000*(crossprod(vecchia_blocks_mult(x = w, vecchia_blocks = vecchia_blocks_, n_time = n_time)) - 
  crossprod(vecchia_blocks_mult(x = w, vecchia_blocks = vecchia_blocks, n_time = n_time)))
}

plot(derivative_mat[lower.tri(derivative_mat)],
tatato)





Sig  = GpGp::exponential_isotropic(c(1, 1, 0), matrix(runif(2000), 1000))
Rho = matrix(1, 1000, 1000)
Rho[seq(500), -seq(500)] = .5
Rho[-seq(500), seq(500)] = .5
Rho_ = matrix(1, 1000, 1000)
Rho_[seq(500), -seq(500)] = .5001
Rho_[-seq(500), seq(500)] = .5001
z = rnorm(1000)

(z %*% solve(Sig*Rho_) %*% z - z %*% solve(Sig*Rho) %*% z)*10000

d_Rho = matrix(0, 1000, 1000)
d_Rho[seq(500), -seq(500)] = 1
d_Rho[-seq(500), seq(500)] = 1

-z %*% solve(Sig*Rho) %*% (d_Rho*Sig) %*% solve(Sig*Rho) %*% z
