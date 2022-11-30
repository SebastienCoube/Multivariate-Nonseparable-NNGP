#guinness_groups = function(NNarray)
#{
#  NNlist = GpGp::group_obs(NNarray, exponent = 3)
#  blocks = split(
#    NNlist$all_inds, 
#    findInterval(seq(length(NNlist$all_inds)), NNlist$last_ind_of_block+1)
#  )
#  resp_in_blocks = split(
#    NNlist$local_resp_inds, 
#    findInterval(seq(length(NNlist$local_resp_inds)), NNlist$last_resp_of_block+1)
#  )
#  return(list(
#    "children"= mapply(function(block, resp)block[resp],  blocks, resp_in_blocks), 
#    "parents" = mapply(function(block, resp)block[-resp], blocks, resp_in_blocks)
#      ))
#}
