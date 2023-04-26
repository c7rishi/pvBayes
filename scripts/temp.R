
library(pvLRT)
tmp <- pvbayes(statin46)
samp_rvar <- stan_output_to_rvar(tmp)
library(posterior)
mean(samp_rvar)
Pr(samp_rvar<1)



class(samp_rvar)

samp_rvar[1,2]



