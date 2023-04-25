# computes cross-product of array along the first and third dimension (typically gives a crossprod per variable)
array_crossprod = function(array_1, array_2 = NULL)
{
  if(is.null(array_2))
  {
    res =  matrix(0, dim(array_1)[2], dim(array_1)[2])
    for(i in seq(ncol(res)))
    {
      for(j in seq(1, i))
      {
        res[i, j] = sum(na.omit(c(array_1[,i,]*array_1[,j,])))
        res[j, i] = res[i, j]
      }
    }
  }
  if(!is.null(array_2))
  {
    res =  matrix(0, dim(array_1)[2], dim(array_2)[2])
    for(i in seq(nrow(res)))
    {
      for(j in seq(ncol(res)))
      {
        res[i, j] = sum(na.omit(c(array_1[,i,]*array_2[,j,])))
      }
    }
  }
  res
}

# computes product of array and matrix following the second dimension of the array and the first dimension of the matrix
# array-matrix equivalent of matrix-vector X %*% beta
# a of size j,k,l,   m of size k,p, a*m of size j,p,l
array_matrix_mult = function(a, m)
{
  res = array(0, dim = c(dim(a)[1], dim(m)[2], dim(a)[3]))
  for(i in seq(nrow(m)))for(j in seq(ncol(m)))res[,j,] = res[,j,] + a[,i,]*m[i,j]
  res
}