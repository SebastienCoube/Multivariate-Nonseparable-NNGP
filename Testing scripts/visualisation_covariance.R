a = .5 # between 0 and 1
b = .9 # between 0 and 1
lambda = 1 # between 0 and 1
r = .0001 # bigger than 0
cc = 10 # bigger than 0
A  = c(1, 1, 1) # between 0 and 1

i=1
j=1

u = seq(0, 1, .01)

plot(u, 
     1/
       (
         (1+cc*u^(2*a))^b-
           A[i]*A[j]*(1+r*u^(2*lambda))^-b
       )
     )
