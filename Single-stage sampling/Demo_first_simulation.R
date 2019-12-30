# This is a demo code to illustrate the first simulation in the paper 'Matrix
# Completion with Covariate Information'
rm(list = ls())
source("SMC_functions.R")  # Load the auxiliary functions, and see Readme.txt for details. 

############################################################## Parameters ##########################
N.popu = 10000  # Population size, the value of $N$ in the paper.
N.ques = 500  # Number of questions, the value of $L$ in the paper.
N.covariates = 20  # Number of covariates, the value of $d$ in the paper.
n.low.rank = 10  # Rank of the response matrix $B_N^*$ in the paper.
informative.ind = FALSE  # The default is non-informative sampling. Change to informative.ind = TRUE for informative sampling.

############################################################## Generation of the finite population #############
set.seed(1234)
beta = matrix(rnorm(N.ques * N.covariates, mean = 0.5, sd = 1), nrow = N.covariates)  # Generate the regression parameter matrix $\beta_N$
beta[beta < quantile(beta, 0.5)] = 0  # Smaller half values of $\b_N^*$ are set to 0.
X.cov = matrix(rnorm(N.popu * N.covariates, mean = 0.5, sd = 1), nrow = N.popu)  # Generate the covariate matrix $X_N$.
BL = matrix(rnorm(N.popu * n.low.rank, mean = 1, sd = 3), nrow = N.popu)  # Generate  $B_L$
BR = matrix(rnorm(N.ques * n.low.rank, mean = 1, sd = 3), nrow = N.ques)  # Generate $B_R$
PX = X.cov %*% solve(t(X.cov) %*% X.cov) %*% t(X.cov)  # Generate the projection matrix to guarantee the orthogonality of $X_N$ and $B_N^*$ on the population level.
Eye = diag(1, dim(X.cov)[1])  # Generate the identity matrix
PXp = Eye - PX
B0 = PXp %*% BL %*% t(BR)
A0 = X.cov %*% beta + B0

############################################################## Generation of missing probability ##############
gamma.i.1 = matrix(rnorm(3 * N.ques, -0.1, 0.1), nrow = 3)  # Generate $\gamma_j$ for $j=1,\ldots,3$.
gamma.i.0 = matrix(rnorm(1 * N.ques, 0.3, 0.1), nrow = 1)  # Generate $\gamma_0$
prob.mis = 1/(1 + exp(-(cbind(1, X.cov[, 1:3]) %*% rbind(gamma.i.0, gamma.i.1))))  # Only the first three components are used for the response model.

############################################################## Generation of inclusion probabilities ############ For non-informative sampling
if (!informative.ind) {
  
  auxiliary.x = apply((X.cov), 1, mean)
  
  size.measure = abs(auxiliary.x) + rexp(N.popu) + 1
  pps.prob = size.measure/sum(size.measure)  # Generate the inclusion probability $\pi_i$.
} 

Y = A0 + matrix(rnorm(N.popu * N.ques, 0, 12), nrow = N.popu)  # Generate the finite population of interest.
target = apply(Y, 2, mean)  # The parameter of interest, that is, the population mean of each question.

if (informative.ind) {
  # For informative sampling
  Y.mean = apply(Y[, 1:7], 1, mean)
  size.measure = ((Y.mean - min(Y.mean)) + 1)
  pps.prob = size.measure/sum(size.measure)  # Generate the inclusion probability $\pi_i$.
}


design = "POI"  # Set the sampling design to be Poisson sampling, and change it to ``PPS'' or ``SRS'' to get the result for probability-proportional-to-size sampling or simple random sampling, respectively.
sample.size = 200  # Set sample size $n=200$, and change the value to 500 get the results for $n=500$.    
cv_ratio_SMC = SMCfit_cv(Y, X.cov, prob.mis, pps.prob, sample.size, design, nfolds = 5, 
    tau1_grid = seq(0, 2, length = 30), tau2_grid = seq(0.9, 0.1, length = 30), alpha_grid = seq(0.992, 
        1, length = 20))  # Choose the tuning parameter by cross-validation.
seed_num = 1
one.iter(seed_num = seed_num, Y = Y, X.cov = X.cov, prob.mis = prob.mis, 
    pps.prob = pps.prob, sample.size = sample.size, design = design, cv_ratio_SMC)
