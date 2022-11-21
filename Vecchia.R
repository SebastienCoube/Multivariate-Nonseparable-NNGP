
vecchia_multi = function(locs, NNarray, var_tag, u, 
                         a, b, cc, delta, lambda, r, 
                         A_vec, nu_vec, a2_vec, rho_vec)
{
  if(nrow(locs)!= length(var_tag))stop("lengths of locs and var_tag do not match")
  
  multiplier_and_variations = get_multiplier_and_variations
  (
    a, b, cc, delta, lambda, r, 
    A_vec, nu_vec, a2_vec, rho_vec, 
    u
  )
  
  effective_range_and_variations = get_effective_range_and_variations
  (
    a, b, cc, delta, lambda, r, 
    A_vec, nu_vec, a2_vec,  
    u
  )
  
  for(i in seq(2, nrow(NNarray)))
  {
    idx = na.omit(NNarray[i,])
    
  }
}

