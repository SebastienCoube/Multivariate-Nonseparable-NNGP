## # generating data 
## n_loc = 500
## n_time = 200
## n_var = 5
## 
## # full data
## y = array(dim = c(n_loc, n_var, n_time))
## dimnames(y) = list(
##   paste("loc", seq(n_loc), sep = "_"), 
##   paste("var", seq(n_var), sep = "_"), 
##   paste("time", seq(n_time), sep = "_")
##   )
## y[] = rnorm(length(y))
## 
## X = array(dim = c(n_loc, 10, n_time))
## dimnames(X) = list(
##   paste("loc", seq(n_loc), sep = "_"), 
##   paste("covariate", seq(dim(X)[2]), sep = "_"), 
##   paste("time", seq(n_time), sep = "_")
##   )
## X[] = rnorm(length(X))
## X[,1,] = 1
## X_noise = X
## X_scale = X
## 
## locs = matrix(runif(2*n_loc), n_loc)


#  symmetric matrix
get_symmetric_mat = function(mat_coordinates)
{
  M = matrix(0,nrow_square_matrix(length(mat_coordinates)),nrow_square_matrix(length(mat_coordinates)))
  M[lower.tri(M, diag = T)] = mat_coordinates
  M = t(M)
  M[lower.tri(M, diag = T)] = mat_coordinates
  M
}
nrow_square_matrix = function(mat_coord_length)(-1+ sqrt(1 + 8*mat_coord_length))/2


# get exponential of symmetric matrix and its derivatives 
get_expmat_and_derivatives = function(mat_coordinates)
{ 
  M_0 = get_symmetric_mat(mat_coordinates)
  eig = eigen(M_0)
  # using directional derivative formula
  G = outer(eig$val, eig$val, function(x, y)(exp(x)-exp(y))/(x-y))
  diag(G)=exp(eig$val)
  G[is.nan(G)]=0
  J = get_lower_tri_idx(ncol(M_0), T)
  list(
    "expmat" = t(eig$vectors)%*%(exp(eig$values)*eig$vectors),
    "expmat_derivative" = lapply(seq(length(mat_coordinates)), function(i)
    {
      #coord = rep(0, length(mat_coordinates)); coord[i]=1; V = get_symmetric_mat(coord) 
      #eig$vectors %*% ((t(eig$vec) %*% V %*% eig$vectors)*G)%*% t(eig$vectors)
      cross  = tcrossprod(eig$vectors[J[i,1],], eig$vectors[J[i,2],])
      #cross+t(cross)-t(eig$vec) %*% V %*% eig$vectors
      eig$vectors %*% ((cross+t(cross))*G)%*% t(eig$vectors)
    })
  )
}

get_tau_precision = function(expmat, y)solve(expmat[!is.na(y), !is.na(y)])

# get the noise precision matrix and its derivatives following NA pattern of y 
get_tau_info = function(mat_coordinates, y)
{
  if(all(is.na(y)))return(NULL)
  if(all(is.na(mat_coordinates)))return(NULL)
  res = list()
  expmat_and_derivatives = get_expmat_and_derivatives(mat_coordinates = mat_coordinates)
  res$tau_precision = get_tau_precision(expmat_and_derivatives$expmat, y)
  res$tau_precision_deriv = lapply(expmat_and_derivatives$expmat_derivative, function(x) res$tau_precision %*% x[!is.na(y), !is.na(y)] %*% res$tau_precision)
  return(res)
}

# get all noise precision matrices and their derivatives following NA pattern of y 
get_noise_info = function(X_noise_list, noise_beta, y_NA_possibilities_match, y_NA_possibilities, y)
{
  mat_coordinates_array = array_matrix_mult(X_noise_list$X, noise_beta)
  # nonstationary case
  if(!X_noise_list$is_only_intercept) expmat_and_derivatives = parallel::mclapply(
    seq(dim(y)[1]), function(loc_index)
    {
      lapply(seq(dim(y)[3]), function(time_index)
      {
        return(get_tau_info(
          mat_coordinates = X_noise_list$X[loc_index,,time_index]%*%noise_beta,
          y = y[loc_index,,time_index]
        ))
      })
    }, mc.cores = parallel::detectCores()-1)
  # stationary case
  if(X_noise_list$is_only_intercept)
  {
    tau_info_cases = lapply(
      split(y_NA_possibilities, row(y_NA_possibilities)), 
      function(x)
      {
        x[x==F]=NA
        get_tau_info(mat_coordinates = rep(1, nrow(noise_beta))%*%noise_beta, y = x)
      }
    )
    expmat_and_derivatives = lapply(seq(dim(y)[1]), function(loc_index)
    {
      lapply(seq(dim(y)[3]), function(time_index)
      {
        tau_info_cases[[y_NA_possibilities_match[loc_index, time_index]]]
      })
    })
  }
  names(expmat_and_derivatives) = paste("loc_", seq(length(expmat_and_derivatives)), sep = "")
  for(x in expmat_and_derivatives)names(x)=paste("time_", seq(length(x)), sep = "")
  return(expmat_and_derivatives)
}


#mat_coordinates = rnorm(15) 
#expmat_and_derivatives = get_expmat_and_derivatives(mat_coordinates)
#y = rnorm(5); y[c(3, 5)] = NA
#
#-.5 * log(det(expmat_and_derivatives$expmat[c(1, 2, 4), c(1, 2, 4)]))-
#  .5 * na.omit(y)%*%solve(expmat_and_derivatives$expmat[c(1, 2, 4), c(1, 2, 4)])%*%na.omit(y)



get_noise_precisions = function(
    noise_info, 
    useful_stuff, 
    Vecchia_approx_DAG, 
    time_begin, 
    time_end
)
{
  # time periods
  k = seq(time_begin, time_end)
  noise_precision_matrices = 
  #  parallel::mc
  lapply(
      #mc.cores = parallel::detectCores()-1,
      X = k, 
      function(i_time)
      {   
        noise_precision_x = 
          unlist(lapply(seq(useful_stuff$n_loc), function(i_loc)
          {
            noise_info[[i_loc]][[i_time]]$tau_precision
          }))
        if(is.null(noise_precision_x))return(
          Matrix::sparseMatrix(
            i = 1,
            j = 1, 
            x = 0, 
            dims = c(useful_stuff$n_field, useful_stuff$n_field)
          )
        )
        Matrix::sparseMatrix(
          i = c(useful_stuff$noise_precision_i[[i_time]]),
          j = c(useful_stuff$noise_precision_j[[i_time]]), 
          x = c(noise_precision_x), 
          dims = c(useful_stuff$n_field, useful_stuff$n_field)
        )
      })
  
}