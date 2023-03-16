target  = .2
s = 0
eps = .001
for(i in seq(100))
{
  p_accepted = 1/(1+exp(- (102 - 103 * s)))
  accepted = rbinom(1,1,p_accepted)
  s = s + accepted * eps * (1-target)
  s = s  - (1-accepted) * eps * target
  print(s)
  print(p_accepted)
}

target  = .2
s = 0
eps = .1
for(i in seq(100))
{
  p_accepted = 1/(1+exp(- (-10 - 1 * s)))
  accepted = rbinom(1,1,p_accepted)
  s = s + accepted * eps * (1-target)
  s = s  - (1-accepted) * eps * target
  print(s)
  print(p_accepted)
}



homesostasy  = function(log_var, accepted, target  = .25, eps = .01)
{
  if(accepted) log_var = log_var + eps * (1-target)
  if(!accepted) log_var = log_var - eps * target
  return(log_var)
}


log_var = 0
p_accepted_record = rep(0, 200)
for(i in seq(200))
{
  p_accepted = 1/(1+exp(- (-30 - 5 * log_var)))
  accepted = rbinom(1,1,p_accepted)
  log_var = homesostasy(log_var, accepted, eps = 10/i)
  print(p_accepted)
  p_accepted_record[i]=p_accepted
}
plot(p_accepted_record)


log_var = 0
p_accepted_record = rep(0, 1000)
for(i in seq(1000))
{
  p_accepted = 1/(1+exp(- (+30 - 5 * log_var)))
  accepted = rbinom(1,1,p_accepted)
  log_var = homesostasy(log_var, accepted, eps = 5/i)
  print(p_accepted)
  p_accepted_record[i]=p_accepted
}
plot(p_accepted_record)






