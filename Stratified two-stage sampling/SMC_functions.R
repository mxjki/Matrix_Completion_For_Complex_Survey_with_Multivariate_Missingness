
########################## Validation error
cv_validate_sim = function(diagD, tau2_grid, tau1_grid_length, svdY, alpha, uvalidatemat, 
                           Xbeta_validate, umat_validate) {
  
  # This function calculates the validation error to tuning parameters.
  
  cv_validate = array(0, c(length(tau2_grid), tau1_grid_length))
  
  for (i in 1:length(tau2_grid)) {
    
    B_validate = as.vector(sqrt(diagD) * SVTE_alpha(svdY$u, svdY$d, svdY$v, tau2_grid[i], 
                                                    alpha))[uvalidatemat != 0]
    
    for (j in 1:tau1_grid_length) {
      
      cv_validate[i, j] = sum((Xbeta_validate[j, ] + B_validate - umat_validate)^2)/sqrt(sum(uvalidatemat != 
                                                                                               0))
      
    }
  }
  
  return(cv_validate)
}



########################## Scale SVT
SVTE_alpha = function(svdu, svdd, svdv, tau_svd_ratio, alpha_ratio) {
  
  # This function performs the singular value soft-thresholding and scaling
  # procedures.
  
  return(svdu %*% (pmax(svdd - svdd[1] * tau_svd_ratio * alpha_ratio, 0) * t(svdv))/(1 + 
                                                                                       2 * svdd[1] * tau_svd_ratio * (1 - alpha_ratio)))
  
}

########################## After obtain the estimator from SMCE, we need transform back by

### Ahat=diagD^{1/2}*Bhat

######################################## Proposed Method
SMCfit_cv = function(Y, X.cov, prob.mis, pps.prob, sample.size,  nfolds = 5, 
                     tau1_grid = seq(0, 1, length = 30), tau2_grid = seq(0.9, 0.1, length = 30), 
                     alpha_grid = seq(0.992, 1, length = 20),sampling.str=sampling.str,data.popu.df = data.popu.df) {
  
  #' This function performs the matrix completion part in the proposed method with the tunning parameter choosed by cross-validation.
  #'
  #' Arguments:
  #'    Y:                The population matrix containing the values of interest.
  #'    X.cov:            The population matrix of covariate.
  #'    prob.mis:         The missing probability for each element of the population.
  #'    pps.prob:         The normalized inclusion probability for probability-proportional-to-size sampling.
  #'    sample.size:      The desired sample size.
  #'    nfolds:           The number of cross-validation folds; with default value 5.
  #'    tau1_grid:        The grid used in search for parameter '$\tau_{1}$'; with default value seq(0,1,length=30).
  #'    tau2_grid:        The grid used in search for parameter '$\tau_{2}$'; with default value seq(0.9,0.1,length=30).
  #'    alpha_grid:       The grid used for alpha; with default value seq(0.992,1,length=20).
  #'
  #' Outputs: 
  #'    A list is returned, with values: Estimated Ahat, Bhat, betahat and rank.
  
  
  N.popu = nrow(Y)
  
  index.s = sample.gen(sample.size = sample.size,sampling.str=sampling.str,data.popu.df = data.popu.df)
  
  sample.full = Y[index.s, ]
  sample.prob.mis = prob.mis[index.s, ]
  x.sample = X.cov[index.s, ]
  x.prob = pps.prob[index.s]
  weight.i.hat = 1/x.prob
  
  sample.mis.result = sample.miss(sample.full, sample.prob.mis = sample.prob.mis)

  diagD = sample.size * x.prob
  
  omega = matrix(as.double(!is.na(sample.mis.result)), nrow = nrow(sample.mis.result))
  
  missing.index = data.frame(index = which(!is.na(sample.mis.result)))
  missing.index$group = sample(x = 1:5, size = nrow(missing.index), replace = TRUE, 
                               prob = rep(1, 5))
  
  Aobs = Y[index.s, ] * omega
  n1 = dim(Aobs)[1]
  n2 = dim(Aobs)[2]
  m = dim(x.sample)[2]
  
  Xtheta = x.sample[, 1:3]
  
  Xadj = x.sample/sqrt(diagD)
  PX = Xadj %*% solve(t(Xadj) %*% Xadj) %*% t(Xadj)
  Eye = diag(1, n1)
  PXp = Eye - PX
  
  iomega = as.numeric(omega != 0)
  iobs = which(iomega == 1)
  imiss = which(iomega == 0)
  
  gammahat = matrix(rep(NA, n1 * n2), nrow = n1)
  for (i in 1:n2) {
    
    thetadata = data.frame(cbind(omega[, i], Xtheta))
    thetaglmout = glm(X1 ~ ., family = binomial(logit), data = thetadata)
    gammahat[, i] = predict(thetaglmout, type = "response")
    
  }
  
  Aobsw = (Aobs/gammahat)/sqrt(diagD)
  
  PXpAobsw = PXp %*% Aobsw
  
  svdPXpAobsw = svd(PXpAobsw)
  
  folds <- cut(sample(seq(1, length(iobs))), breaks = nfolds, labels = FALSE)
  iomegatrain = matrix(0, n1 * n2, nfolds)
  iomegavalidate = matrix(0, n1 * n2, nfolds)
  omegatrain = array(NA, c(n1, n2, nfolds))
  omegavalidate = array(NA, c(n1, n2, nfolds))
  
  gammatrainhat = array(NA, c(n1, n2, nfolds))
  Aobswtrain = array(NA, c(n1, n2, nfolds))
  PXpAobswtrain = array(NA, c(n1, n2, nfolds))
  svdPXpAobswtrain = list()
  A_validate = list()
  Xbetavalidate = list()
  
  for (fold_num in 1:nfolds) {
    
    validateIndexes <- which(folds == fold_num, arr.ind = TRUE)
    ivalidate <- iobs[validateIndexes]
    itrain <- iobs[-validateIndexes]
    
    iomegatrain[itrain, fold_num] = 1
    omegatrain[, , fold_num] = matrix(iomegatrain[, fold_num], n1, n2)
    
    iomegavalidate[ivalidate, fold_num] = 1
    omegavalidate[, , fold_num] = matrix(iomegavalidate[, fold_num], n1, n2)
    
    for (i in 1:n2) {
      
      thetadata = data.frame(cbind(omegatrain[, i, fold_num], Xtheta))
      thetaglmout = glm(X1 ~ ., family = binomial(logit), data = thetadata)
      gammatrainhat[, i, fold_num] = predict(thetaglmout, type = "response")
      
    }
    
    Aobswtrain[, , fold_num] = omegatrain[, , fold_num] * Aobs/gammatrainhat[, 
                                                                             , fold_num]
    PXpAobswtrain[, , fold_num] = PXp %*% Aobswtrain[, , fold_num]
    svdPXpAobswtrain = append(svdPXpAobswtrain, list(svd(PXpAobswtrain[, , fold_num])))
    
    A_validate = append(A_validate, list(as.vector(Aobs)[iomegavalidate[, fold_num] != 
                                                           0]))
    
    Xbetavalidate = append(Xbetavalidate, list(matrix(nrow = length(tau1_grid), 
                                                      ncol = sum(iomegavalidate[, fold_num]))))
  }
  
  for (j in 1:length(tau1_grid)) {
    
    PX = Xadj %*% solve(t(Xadj) %*% Xadj + svd(t(Xadj) %*% Xadj)$d[1] * tau1_grid[j] * 
                          diag(1, m)) %*% t(Xadj)
    
    for (fold_num in 1:nfolds) {
      Xbetavalidate[[fold_num]][j, ] = as.vector(PX %*% Aobswtrain[, , fold_num])[iomegavalidate[, 
                                                                                                 fold_num] != 0]
    }
  }
  
  cv = array(0, c(nfolds, length(tau2_grid), length(tau1_grid), length(alpha_grid)))
  for (fold_num in 1:nfolds) {
    for (alpha_num in 1:length(alpha_grid)) {
      
      cv[fold_num, , , alpha_num] = cv_validate_sim(diagD, tau2_grid, length(tau1_grid), 
                                                    svdPXpAobswtrain[[fold_num]], alpha_grid[alpha_num], iomegavalidate[, 
                                                                                                                        fold_num], Xbetavalidate[[fold_num]], A_validate[[fold_num]])
      
    }
  }
  
  cv_grid = which(apply(cv, 2:4, sum) == min(apply(cv, 2:4, sum)), arr.ind = T)[1, 
                                                                                ]
  
  return(list(cv = cv, tau_beta_ratio = tau1_grid[cv_grid[2]], tau_svd_ratio = tau2_grid[cv_grid[1]], 
              alpha_ratio = alpha_grid[cv_grid[3]]))
}


SMCfit = function(Aobs, X, diagD, tau_beta_ratio, tau_svd_ratio, alpha_ratio) {
  
  #' This function performs the matrix completion part in proposed method with fix tunning parameters.
  #'
  #' Arguments:
  #'    Aobs:             The observation sample matrix.
  #'    X:                The sample covariate matrix.
  #'    diagD:            The input diagonal matrix with inclusion probabilities pi_{i}.
  #'    tau_beta_ratio:   The input tunning parameter 'tau_{1}'.
  #'    tau_svd_ratio:    The input tunning parameter 'tau_{2}'.
  #'    alpha_ratio:      The tunning parameter 'alpha'.
  #'
  #' Outputs: 
  #'    A list is returned, with values: Estimated Ahat, Bhat, betahat and rank.
  
  n1 = dim(Aobs)[1]
  n2 = dim(Aobs)[2]
  m = dim(X)[2]
  
  iomega = as.numeric(Aobs != 0)
  omega = matrix(iomega, n1, n2)
  
  Xtheta = X[, 1:3]
  
  Xadj = X/sqrt(diagD)
  PX = Xadj %*% solve(t(Xadj) %*% Xadj) %*% t(Xadj)
  Eye = diag(1, n1)
  PXp = Eye - PX
  
  gammahat = matrix(rep(NA, n1 * n2), nrow = n1)
  for (i in 1:n2) {
    
    thetadata = data.frame(cbind(omega[, i], Xtheta))
    thetaglmout = glm(X1 ~ ., family = binomial(logit), data = thetadata)
    gammahat[, i] = predict(thetaglmout, type = "response")
    
  }
  
  Aobsw = (Aobs/gammahat)/sqrt(diagD)
  
  PXpAobsw = PXp %*% Aobsw
  
  ### The following part is different from the original one.
  betahat = solve(t(Xadj) %*% Xadj + svd(t(Xadj) %*% Xadj)$d[1] * tau_beta_ratio * 
                    diag(1, m)) %*% t(Xadj) %*% Aobsw
  Xbetahat = X %*% betahat
  
  svdPXpAobsw = svd(PXpAobsw)
  
  Bhat = SVTE_alpha(svdPXpAobsw$u, svdPXpAobsw$d, svdPXpAobsw$v, tau_svd_ratio, 
                    alpha_ratio)
  
  rankB = sum(pmax(svdPXpAobsw$d - svdPXpAobsw$d[1] * tau_svd_ratio * alpha_ratio, 
                   0) > 0)
  
  ###################### 
  Ahat = Xbetahat + Bhat
  
  return(list(Ahat = Ahat, Bhat = Bhat, betahat = betahat, rank = rankB + m))
}


######################################## Data generation ##########################

# auxiliary.x = auxiliary.c We only consider MCAR
sample.gen = function(sample.size, sampling.str,data.popu.df) {
  #' This function is used to generate a sample.   
  #' 
  #' Arguments:
  #'    sample.size:  The desired sample size within each cluster.
  #'    sampling.str: The stratum and cluster information for each element.
  #'    data.popu.df: This data frame contains the number of elements within each cluster.
  #'                   
  #' Outputs:
  #'    The index of the sample.    
   
   clust.select = data.frame(
      strata = rep(1:max(data.popu.df$strata),each = 2),
      clst = NA
    )
    for (i in 1:max(data.popu.df$strata))
    {
      clus.popu = subset(data.popu.df,strata == i)
      clust.select$clst[clust.select$strata == i] = sample(x = 1:nrow(clus.popu),2,prob = clus.popu$c.size,replace = TRUE)
    }
    index.s = rep(NA,max(data.popu.df$strata)*2* sample.size)
    
    index.j.i = 1
    for (j in 1:nrow(clust.select)){
      popu.j = subset(sampling.str,strata == clust.select$strata[j] & cluster == clust.select$clst[j] )
      index.s[ index.j.i:  (index.j.i + sample.size-1)] = sample(popu.j$index, sample.size,replace = FALSE)
      index.j.i = index.j.i +  sample.size
    }
    return(sort(index.s))
}

sample.miss = function(sample.full,sample.prob.mis) {
  #' This function is used to generate the missingness for the sample.
  #'    
  #' Arguments:
  #'    sample.full:     The sample without missing. 
  #'    sample.prob.mis: The missing probability for each element of the sample.
  #'    
  #' Outputs:
  #'    The sample with missing value, denoted as 'NA'.
    
    sample.vec = c(sample.full)
    miss.ind = rbinom(n = length(sample.vec), size = 1, prob = c(sample.prob.mis))
    
    sample.vec[which(miss.ind == 0)] = NA
    sample.mis = matrix(sample.vec, nrow = nrow(sample.full))
    return(sample.mis)
}


est.double.r = function(fitted.sample, ori.sample, est.p,  pps.prob, N.popu, sample.size) {
  
  #' This function is used to apply the doubly robust estimator.   
  #' 
  #' Arguments:
  #'    fitted.sample: The fitted sample using a specific model.
  #'    ori.sample:    The original sample with missing data denoted as 'NA'.  
  #'    est.p:         The sample without missing for a specific column. 
  #'    pps.prob:      The normalized inclusion probability with respect to each stratum.
  #'    N.popu:        The population size
  #'    sample.size:   The desired sample size.
  #'    
  #' Outputs:
  #'    Estimated population mean using a doubly robust estimator. 
  
  omega.i = as.double(ori.sample!=0)
  element.part = fitted.sample
  element.part[omega.i == 1] = element.part[omega.i == 1] + ((ori.sample - fitted.sample)/est.p)[omega.i == 1]
  
  est.result = sum(element.part/pps.prob)/N.popu/2/sample.size
  
  
  return(est.result)
}

est.double.r.v = function(fitted.sample, ori.sample, est.p, diagD.ii, N.popu, sample.size) {

  #' This function is used to apply the doubly robust estimator for the variance estimation.   
  #' 
  #' Arguments:
  #'    fitted.sample: The fitted sample using a specific model.
  #'    ori.sample:    The original sample with missing data denoted as 'NA'.  
  #'    est.p:         The sample without missing for a specific column. 
  #'    diagD.ii:      The replicated inclusion probability for the modified repeated replication method.
  #'    N.popu:        The population size
  #'    sample.size:   The desired sample size.
  #'    
  #' Outputs:
  #'    Estimated population mean using a doubly robust estimator.
  
  omega.i = as.double(ori.sample!=0)
  element.part = fitted.sample
  element.part[omega.i == 1] = element.part[omega.i == 1] + ((ori.sample - fitted.sample)/est.p)[omega.i == 1]
  
  est.result = sum(element.part/diagD.ii)/N.popu
  
  return(est.result)
}
one.iter = function(seed_num, Y, X.cov, prob.mis, pps.prob, sample.size = sample.size, cv_ratio_SMC,
                    sampling.str=sampling.str,data.popu.df = data.popu.df,hadm = hadm) {
  #' This function is used to estimate the population mean by five different methods and the variance of the doubly robust
  #' estimator using the proposed method.
  #'    
  #' Arguments:
  #'    seed_num:         The value used for the random seed.
  #'    Y:                The population matrix containing the values of interest.
  #'    X.cov:            The population matrix of covariate.
  #'    prob.mis:         The missing probability of elements in Y.
  #'    pps.prob:         The normalized inclusion probability with respect to each stratum.
  #'    sample.size:      The desired sample size.
  #'    design:           Sampling design used to get the sample.  
  #'    cv_ratio_SMC:   Tuning parameter for the proposed matrix completion method.
  #'
  #'  Outputs:
  #'    Error of population mean estimator using five methods and the plug-in variance estimator for the doubly robust estimator using the proposed method. 
    
    set.seed(1564 + seed_num * 3)
    N.popu = nrow(Y)
    target = apply(Y, 2, mean)

    index.s = sample.gen(sample.size = sample.size, sampling.str=sampling.str,data.popu.df = data.popu.df)
    sample.full = Y[index.s, ]
    sample.prob.mis = prob.mis[index.s, ]
    x.sample = X.cov[index.s, ]
    x.prob = pps.prob[index.s]
    sample.mis.result = sample.miss(sample.full, sample.prob.mis = sample.prob.mis)
    diagD = sample.size * x.prob * 2
    
    omega = matrix(as.double(!is.na(sample.mis.result)), nrow = nrow(sample.mis.result))
    Phat = sample.mis.result * 0 + NA
    for (ques.i in 1:ncol(sample.mis.result)) {
      
      y.i.hat = c(omega[, ques.i])
      weight.i.hat = 1/x.prob
      
      est.i.hat = glm(y.i.hat ~ x.sample[, 1:3], weights = weight.i.hat/2000, family = binomial(link = "logit"))$fitted.values
      Phat[, ques.i] = est.i.hat
    }
      
   ind.na.i = is.na(sample.mis.result)
   sample.mis.result[ind.na.i] = 0
   sample.mis.SMC = SMCfit(Aobs = sample.mis.result, X = x.sample, diagD, cv_ratio_SMC$tau_beta_ratio, 
                             cv_ratio_SMC$tau_svd_ratio, cv_ratio_SMC$alpha_ratio)
   fitted.beta.sample.SMC = sample.mis.SMC$beta
   fitted.B.sample.SMC = sqrt(diagD) * sample.mis.SMC$Bhat
   fitted.sample.SMC = x.sample %*% sample.mis.SMC$beta + sqrt(diagD) * sample.mis.SMC$Bhat
   
   est.result <- var.result <- matrix(NA, nrow = 1, ncol = ncol(Y))
   for (i in 1:(ncol(Y))) {
     est.result[i] = est.double.r(fitted.sample  = fitted.sample.SMC[, i], ori.sample  = sample.mis.result[, i], est.p = Phat[, i], pps.prob = x.prob, N.popu = N.popu, sample.size = sample.size) - target[i]
   }
   
   boo.rep = matrix(NA,nrow = nrow(hadm), ncol =ncol(Y))
   for (ii in 1:nrow(hadm)){
     had.ii = hadm[ii,]
     fac.ii = rep(had.ii,each = 2)
     fac.ii[2*(1:length(had.ii))] = -had.ii
     
     ## We consider the BRR(1/2)
     x.prob.ii = x.prob 
     diagD.ii = ((sample.size * x.prob.ii * 2)^(-1)* rep(1 + 1/2*fac.ii,each = sample.size))^(-1) 
     
     Phat = sample.mis.result * 0 + NA
     for (ques.i in 1:ncol(sample.mis.result)) {
       
       y.i.hat = c(omega[, ques.i])
       weight.i.hat = 1/x.prob.ii*rep(1 + 1/2*fac.ii,each = sample.size)
       
       est.i.hat = glm(y.i.hat ~ x.sample[, 1:3], weights = weight.i.hat/2000, family = binomial(link = "logit"))$fitted.values
       Phat[, ques.i] = est.i.hat
     }
     
     ind.na.i = is.na(sample.mis.result)
     sample.mis.result[ind.na.i] = 0
     
     sample.mis.SMC = SMCfit(Aobs = sample.mis.result, X = x.sample, diagD.ii, cv_ratio_SMC$tau_beta_ratio, 
                               cv_ratio_SMC$tau_svd_ratio, cv_ratio_SMC$alpha_ratio)
     fitted.beta.sample.SMC = sample.mis.SMC$beta
     fitted.B.sample.SMC = sqrt(diagD.ii) * sample.mis.SMC$Bhat
     fitted.sample.SMC = x.sample %*% sample.mis.SMC$beta + sqrt(diagD.ii) * sample.mis.SMC$Bh
     
     for (i in 1:(ncol(Y))) {
       boo.rep[ii,i] = est.double.r.v(fitted.sample = fitted.sample.SMC[, i], ori.sample = sample.mis.result[, i], est.p = Phat[, i],  diagD.ii=diagD.ii,N.popu = N.popu, sample.size = sample.size) - target[i]
     }
     
   }
   
   var.result = 4*apply((t(boo.rep) - est.result)^2,1,mean)
   
   return(list(est.proposed = est.result,
               var.proposed = var.result))
    
    
}
