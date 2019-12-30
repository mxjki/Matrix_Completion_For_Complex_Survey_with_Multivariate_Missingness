# This is a demo code to illustrate the second simulation in the paper 'Matrix
# Completion with Covariate Information'
rm(list = ls())
library(plyr)
source("SMC_functions.R")  # Load the auxiliary functions, and see Readme.txt for details. 
hadm = as.matrix(read.csv("Hadamard_15.csv")) # This matrix is used for the modified repeated replication method for variance estimation.

############################################################## Parameters ##########################
H = 15  #Number of strata, the value $H$ in the paper.
N.covariates = 20 # Number of covariates, the value of $d$ in the paper. 
N.popu.basic = 40  # Basic stratum size, and it is used for generating $N_{hi}$. .
N.cluster.basic = 20 # This is the basic cluster size used for generating $N_h$.
N.ques = 500  # Number of questions, the value of $L$ in the paper.
n.low.rank = 10  # Rank of the response matrix $B_N^*$ in the paper.
informative.ind = FALSE  # The default is non-informative sampling. Change to informative.ind = TRUE for informative sampling.

############################################################## Generation of the finite population #############
set.seed(1234)
N.clusters = rpois(H, lambda = N.cluster.basic/2) + N.cluster.basic # Generation of stratum size, the value of $N_h$ in the paper.
zeta.h.single = rexp(H,rate = 3)
zeta.h = rep(zeta.h.single, N.clusters) # The stratum effect, the value of $\zeta_h$ in the paper.  
xi.hi = rexp(sum(N.clusters),rate = 3) # The cluster effect, the value of $\xi_{hi}$ in the paper.
N.cluster.popu = rpois(sum(N.clusters),lambda = N.popu.basic/4*(xi.hi + zeta.h)) + N.popu.basic # The cluster size, the value $N_{hi}$ in the paper.
N.popu = sum(N.cluster.popu) # Population size, the value $N$ in the paper.

data.popu.df = NULL # This data frame contains the cluster size for each stratum, and it is used to generate the finite population.
effect.df = NULL
pop.index = 1
for (i in 1:H){
  data.popu.df = rbind(data.popu.df,
                       data.frame(strata = i,
                                  cluster = 1:N.clusters[i],
                                  c.size =  N.cluster.popu[pop.index:(pop.index+N.clusters[i]-1)]
                       ))
  pop.index = pop.index + N.clusters[i]
}  
sampling.str = data.popu.df[rep(1:nrow(data.popu.df),data.popu.df$c.size),] # This one is used to draw the sample based on a stratified two-stage cluster sampling
sampling.str$index = 1:N.popu
Num.ele.stratum = ddply(data.popu.df,"strata",summarize,size = sum(c.size))$size # This data frame contains number of elements in each stratum.
pps.prob = rep(1/Num.ele.stratum,Num.ele.stratum) # The normalized selection probability with respect to each stratum. 

effects.s.c = rep(zeta.h.single,Num.ele.stratum)/2 + rep(xi.hi,N.cluster.popu)/2  # This one is used to generate the covariate $x_{hi}$.
effects.s.c.2 = rep(zeta.h.single,Num.ele.stratum)/3 + rep(xi.hi,N.cluster.popu)/3 # This one is used for generate the $B_L$.
beta = matrix(rnorm(N.ques * N.covariates, mean = 0.5, sd = 1), nrow = N.covariates) # Generate the regression parameter matrix $\beta_N$
X.cov = matrix(rnorm(N.popu * N.covariates, mean = 0.5, sd = 1), nrow = N.popu) + matrix(rep(effects.s.c,N.covariates),ncol = N.covariates) # Generate the covariate matrix $X_N$.
PX = X.cov %*% solve(t(X.cov) %*% X.cov) %*% t(X.cov) # Generate the projection matrix to guarantee the orthogonality of $X_N$ and $B_N^*$ on the population level.
Eye = diag(1, dim(X.cov)[1])
PXp = Eye - PX

BL = matrix(rnorm(N.popu * n.low.rank, mean =1, sd = 3), nrow = N.popu) + matrix(rep(effects.s.c.2,n.low.rank),ncol = n.low.rank)  # Generate  $B_L$
BR = matrix(rnorm(N.ques * n.low.rank, mean = 1, sd = 3), nrow = N.ques)  # Generate $B_R$
B0 = PXp %*% BL %*% t(BR)
A0 = X.cov %*% beta + B0

############################################################## Generation of missing probability ##############
prob.mis = A0 * 0
gamma.i.1 = matrix(rnorm(3 * N.ques, -0.1, 0.1), nrow = 3) 
gamma.i.0 = matrix(rnorm(1 * N.ques, 0.3, 0.1), nrow = 1) 
prob.mis = 1/(1 + exp(-(cbind(1, X.cov[, 1:3]) %*% rbind(gamma.i.0, gamma.i.1))))

    
Y = A0 + matrix(rnorm(N.popu * N.ques, 0, 12), nrow = N.popu)
target = apply(Y, 2, mean) 
   
sample.size = 10 # Sample Size $n_c$ within each cluster, and change the value to 20 get the results for $n_c=20$.
cv_ratio_SMC = SMCfit_cv(Y, X.cov, prob.mis, pps.prob, sample.size, nfolds = 5, 
                         tau1_grid = seq(0, 2, length = 30), tau2_grid = seq(0.9, 0.1, length = 30), 
                         alpha_grid = seq(0.992, 1, length = 20),sampling.str=sampling.str,
                         data.popu.df = data.popu.df)

i=1
one.iter(seed_num = i,  Y = Y, X.cov = X.cov,  prob.mis = prob.mis, pps.prob = pps.prob, sample.size = sample.size,  cv_ratio_SMC = cv_ratio_SMC,
         sampling.str=sampling.str,data.popu.df = data.popu.df,hadm = hadm)
    