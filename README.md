# Multivariate-Nonseparable-NNGP

An attempt at doing [Nearest Neighbor gaussian Processes](https://arxiv.org/pdf/1406.7343.pdf) with nonseparable space-time [multivariate covariance functions](https://hal.archives-ouvertes.fr/hal-03564931/). The aim is to model air contamination data in the Basque country, it's what I am supposed to do at [BCAM](http://www.bcamath.org/en/) when I'm not eating pintxos de tortilla. 
My approach is to cut corners and save costs in the Vecchia approximation by assuming measure sites that are fixed in space and measurements that are regular in time (there can be missing measurements though). The model will accomodate nonstationarity in the marginal variance of the NNGP and the variance of the Gaussian noise. 
