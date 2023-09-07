# import code for BSIM
source("/gpfs/home/dw2625/r/BSIM/bsim_code.R")
source("/gpfs/home/dw2625/r/BSIM/brms_single_bayes_02152023.R")
#source("./bsim_code.R")
#source("./brms_single_bayes_02152023.R")

library(splines)
library(movMF) 
library(mgcv) 
library(simml)
library(cmdstanr)
library(posterior)
library(data.table)
library(dplyr)
library(slurmR)
library(stableGR)  # to compute the Gelman-Rubin diagnostic statistic

###Generate data
generate.data2 <- function(n = 200, p = 10, #family = "gaussian", 
                           correlationX = 0.2, 
                           sigmaX = 1, #sigma = 0.4, s = 2, delta = 1, 
                           g.choice = "nonlinear", # "linear" 
                           m.choice = "linear",  # "nonlinear"
                           pi.1 = 0.5)
{
  
  true.beta <- c(c(1, 0.5, 0.25, 0.125), rep(0, p - 4))
  true.beta <- true.beta/sqrt(sum(true.beta^2))#sum(true.beta^2)=1
  
  eta.hold <- c(1, 2, 3, 4)
  eta.hold <- eta.hold/sqrt(sum(eta.hold^2))
  true.eta <- c(rep(0, p - 4), eta.hold)
  
  if(g.choice=="nonlinear"){
    g <- function(u) exp(-(u - 0.5)^2) - 0.6 
  }else{
    g <- function(u) 0.3 * u
  }
  
  if(m.choice=="nonlinear"){
    m <- function(u) 0.5* sin(u * 0.5 * pi)
  }else{ # "linear" 
    m <- function(u) 0.125* u * pi
  }
  
  A <- drop(rbinom(n, 1, pi.1) + 1)#A:treatment assignment: 1 or 2
  Psix <- sigmaX * (diag(1 - correlationX, nrow = p, ncol = p) + 
                      matrix(correlationX, nrow = p, ncol = p))#covariance matrix
  ePsix <- eigen(Psix)
  X <- sapply(1:p, function(x) rnorm(n)) %*% diag(sqrt(ePsix$values)) %*% t(ePsix$vectors)
  
  main.effect <- m(drop(X %*% true.eta))
  g.effect <- g(drop(X %*% true.beta))
  interaction.effect <- 2 * (as.numeric(A) + pi.1 - 2) * g.effect
  
  # potential expected outcome 
  mu.1 <- 1/(1 + exp(-(main.effect+  2 * (pi.1 - 1) *g.effect)))#A=1;trt=-0.5
  mu.2 <- 1/(1 + exp(-(main.effect+  2 * pi.1 *g.effect)))#A=2;trt=0.5
  
  # canonical parameter
  theta <- main.effect + interaction.effect
  SNR <- var(interaction.effect)/var(main.effect)
  SNR
  
  mu <- 1/(1 + exp(-theta))#P(Y=1|theta)
  y <- rbinom(n, size = 1, prob = mu)
  
  optTr <- ifelse(mu.2 > mu.1, 1, 2)  # smaller outcome desirable 
  value.opt <- mean((optTr==1)*mu.1 + (optTr==2)*mu.2)
  value.opt ##E(P(Y=1|theta)) using optimal trt
  
  list(y = y, A = A, X = X, SNR = SNR, true.beta = true.beta, 
       true.eta = true.eta, m.choice=m.choice, #delta = delta, s = s, 
       mu.1 = mu.1, mu.2 = mu.2, optTr = optTr, value.opt = value.opt)
}




oneExperimentFn <- function(iter,
                            n = 500, # n=500, 1000, 2000 
                            p = 10,  # p=5, 10
                            g.choice = "linear", # "nonlinear"
                            m.choice = "linear", # "nonlinear" 
                            correlationX = 0.2, 
                            sigmaX = 1,
                            pi.1 = 0.5, 
                            beta.ini = c(1,rep(0,p-1)), # initial value for the frequentist estimate
                            n.test =10000,
                            n.samples= 2000, 
                            n.burning= 2000,
                            lambda.prior=100,
                            lambda.proposal=1000)
{
  
  # generate training data
  data <- generate.data2(n=n, p=p, correlationX=correlationX, sigmaX=sigmaX, g.choice=g.choice,m.choice=m.choice, pi.1=pi.1)
  data$SNR  # the ratio of interactions("signal") vs. main effects("noise")
  A <- data$A
  y <- data$y
  X <- data$X
  Xm <- X
  data$true.beta
  data$true.eta
  # generate test data
  test.data <- generate.data2(n=n.test, p=p, correlationX=correlationX,sigmaX=sigmaX, g.choice=g.choice,m.choice=m.choice, pi.1=pi.1)
  #test.data$y
  xnew <- test.data$X
  xmnew <- xnew   # main effect term
  newA <- test.data$A
  k = 4 + floor(n^{1/5.5})
    
  ## 1) bsim model  
  bsim.obj <-  bsim.fn(y,A,X,Xm, family = "binomial", n.samples=n.samples, 
                       n.burning = n.burning,  
                       lambda.proposal=lambda.proposal,  
                       lambda.prior=lambda.prior, 
                       beta.ini=beta.ini, 
                       k = k)
  bsim.obj$beta.ini 
  #plot(bsim.obj$gamma.sample[3,])
  #plot(bsim.obj$m.sample[3,])
  mean(bsim.obj$acceptsample) 
  apply(bsim.obj$beta.sample, 1, mean)#single index coefficients' posterior mean
  
  # compute the Gelman-Rubin diagnostic statistic
  # GR.beta <- stable.GR(list(t(bsim.obj$beta.sample)))  
  # GR.gamma <- stable.GR(list(t(bsim.obj$gamma.sample)) )  
  # GR.m <- stable.GR(list(t(bsim.obj$m.sample)) )  
  
  # BSIM prediction (apply to the test data)
  predicted <- pred.fn(bsim.obj, xnew=xnew, xmnew=xmnew, newA=newA, thresh=0)
  
  # BSIM estimated optimal treatment decision
  bsim.optim.trt <- sapply(predicted$tbi, function(x) ifelse(x > 0.5, 2, 1)) #if Pr(contrast<0) greater than 0.5, recommend CCP
  
  ## BSIM "value" (expected outcome under ITR)
  value.tmp <- ifelse(bsim.optim.trt==1, test.data$mu.1, test.data$mu.2)
  bsim.value <- mean(value.tmp) #E(P(Y=1|BSIM estimated optimal trt))
  bsim.value 
  value.all.get.1 <- mean(test.data$mu.1) #E(P(Y=1|all get A=1))
  value.all.get.2 <- mean(test.data$mu.2) #E(P(Y=1|all get A=2))
  
  
  # True optimal treatment decision
  true.optim.trt <- test.data$optTr  
  # Percentage of patients' true optimal treatment is CCP 
  
  # The percentage of making correct treatment decision (PCD) using BSIM
  bsim.tab <- table(bsim.optim.trt, true.optim.trt) 
  bsim.accuracy <- sum(diag(bsim.tab))/sum(bsim.tab)
  bsim.accuracy 
  
  
  ## "prediction error" (out-sample deviance)
  # predicted$eta.distr (n.test*n.samples): a canonical parameter matrix with a row for each observation(test data point) and a column for each posterior sample
  ##ll:each row: binomial log-likelihood (for each observation in the test data) use exponential family distribution function
  ll <-  test.data$y *predicted$eta.distr - log(1 + exp(predicted$eta.distr))  # n.test*n.samples;binomial log-likelihood (for the test data) use exponenetial family function
  n.samples <- ncol(ll) 
  f <- function(i){     # this computes log of the mean of exponentiated values for log likelihood for each oberservations in the test data(for numerical stability)
    xmax <- max(ll[i,])
    xsum <- sum(exp(ll[i, ] - xmax)) #ensure that the largest positive exponentiated term is exp(0)=1.
    xmax + log(xsum) - log(n.samples)#mean of likelihood(ie.e, mean of each row of exp(ll)) for observation
  }
  bsim.deviance <- -2*mean(sapply(1:nrow(ll), f)) 
  bsim.deviance
  
  
  
  ### 2) Single Bayesian model
  sgbayes.coefs <- bayes_single_model(X = X,p=p, y = y, A = A) 
  
  #GR.sgbayes <- stable.GR(as.matrix(sgbayes.coefs[,1:(ncol(sgbayes.coefs)-1)]))  
  # Divergent transitions from this Bayesian mdoel
  div <- unique(sgbayes.coefs$div_single)
  sng.bayes.predicted <- pred_sng_bayes(xnew=xnew, xmnew=xmnew, newA=newA, A.mean=bsim.obj$A.mean, coefs = sgbayes.coefs, thresh = 0)
  # Estimated optimal treatment decision
  sng.bayes.optim.trt <- sapply(sng.bayes.predicted$tbi, function(x) ifelse(x > 0.5, 2, 1)) #if Pr(contrast<0) greater than 0.5, recommend CCP
  
  # The percentage of making correct treatment decision using BSIM
  sng.bayes.tab <- table(sng.bayes.optim.trt, true.optim.trt)
  sng.bayes.accuracy <- sum(diag(sng.bayes.tab))/sum(sng.bayes.tab)
  
  ## Single model bayes "value" (expected outcome under ITR)
  sng.bayes.value.tmp <- ifelse(sng.bayes.optim.trt==1, test.data$mu.1, test.data$mu.2)
  sng.bayes.value <- mean(sng.bayes.value.tmp) #E(P(Y=1|single basyes model estimated optimal trt))
  sng.bayes.value 
  
  ## "prediction error" (out-sample deviance)
  # predicted$eta.distr (n.test*n.samples): a canonical paramter matrix with a row for each observation(test data point) and a column for each posterior sample
  ##ll:each row: binomial log-likelihood (for each observation in the test data) use exponential family distributionfunction
  sng.bayes.ll <-  test.data$y *sng.bayes.predicted$eta.distr - log(1 + exp(sng.bayes.predicted$eta.distr))  # n.test*n.samples;binomial log-likelihood (for the test data) use exponenetial family function
  n.samples <- ncol(sng.bayes.ll) 
  
  sng.bayes.f <- function(i){     # this computes log of the mean of exponentiated values for log likelihood for each oberservations in the test data(for numerical stability)
    xmax <- max(sng.bayes.ll[i,])
    xsum <- sum(exp(sng.bayes.ll[i, ] - xmax)) #ensure that the largest positive exponentiated term is exp(0)=1.
    xmax + log(xsum) - log(n.samples)#mean of likelihood(ie.e, mean of each row of exp(ll)) for observation
  }
  
  sng.bayes.deviance <- -2*mean(sapply(1:nrow(sng.bayes.ll), sng.bayes.f)) 
  sng.bayes.deviance
  
  
  data.table(iter=iter,n=n, p=p, g.choice=g.choice,m.choice=m.choice, 
             accept.rate= mean(bsim.obj$acceptsample),
             bsim.value= bsim.value, 
             opt.value = test.data$value.opt,
             diff.value= bsim.value- test.data$value.opt,
             value.all.get.1=value.all.get.1, 
             value.all.get.2=value.all.get.2,
             bsim.deviance=bsim.deviance, 
             bsim.accuracy=bsim.accuracy,
             sng.bayes.accuracy=sng.bayes.accuracy,
             sng.bayes.value =sng.bayes.value,
             sng.bayes.deviance=sng.bayes.deviance,
             sng.bayes.diff.value= sng.bayes.value - test.data$value.opt,
             div=div
             #,
             # GR.beta.psrf.1.03 = mean(GR.beta$psrf<1.03),
             # GR.beta.psrf.1.05 = mean(GR.beta$psrf<1.05),
             # GR.beta.psrf.1.07 = mean(GR.beta$psrf<1.07),
             # GR.gamma.psrf.1.03 =mean(GR.gamma$psrf<1.03),
             # GR.gamma.psrf.1.05 =mean(GR.gamma$psrf<1.05),
             # GR.gamma.psrf.1.07 =mean(GR.gamma$psrf<1.07),
             # GR.m.psrf.1.03 =mean(GR.m$psrf<1.03),
             # GR.m.psrf.1.05 =mean(GR.m$psrf<1.05),
             # GR.m.psrf.1.07 =mean(GR.m$psrf<1.07)
             # #GR.sgbayes.psrf.1.03 = mean(GR.sgbayes$psrf<1.03),
             # GR.sgbayes.psrf.1.05 = mean(GR.sgbayes$psrf<1.05),
             # GR.sgbayes.psrf.1.07 = mean(GR.sgbayes$psrf<1.07)
  )
}




n.test = 10000  # testing set sample size
n.samples=2000
n.burning = 2000
lambda.prior = 100
lambda.proposal =  1000

# test-run with 2 replications for one particular scenario 
# res <- rbindlist(lapply(1:2, function(a) oneExperimentFn(a,n=500, p=5, m.choice="linear", g.choice="nonlinear",
#                                                          lambda.prior=lambda.prior, lambda.proposal=lambda.proposal,
#                                                          n.test=n.test, n.samples=n.samples, n.burning = n.burning)))
# res

scenarios <- expand.grid(n = c(500, 1000, 2000),
                         p = 5,
                         g.choice=c("linear", "nonlinear"),
                         m.choice=c("linear"))
scenarios
results.aggregated2 <- vector("list", length= nrow(scenarios))
n.rep = 100 

set.seed(2022) 
for(r in  1:nrow(scenarios))
{
  print(r)
  n    = scenarios[r, 1]  # training set sample size
  p    = scenarios[r, 2]
  g.choice  = scenarios[r, 3]
  m.choice  = scenarios[r, 4]  
  
  # res <- replicate(1, oneExperimentFn(iter=1,n=n, p=p, g.choice=g.choice, m.choice=m.choice,
  #                                     n.test= n.test,
  #                                     n.samples=n.samples, n.burning = n.burning,
  #                                     lambda.prior=lambda.prior,
  #                                     lambda.proposal=lambda.proposal
  #                                     ))
  # res
  
  sjob <- Slurm_lapply(1:n.rep,
                       FUN=oneExperimentFn,
                       n=n, p=p, g.choice=g.choice, m.choice=m.choice, 
                       n.test= n.test, 
                       n.samples=n.samples, n.burning = n.burning,
                       lambda.prior=lambda.prior,
                       lambda.proposal=lambda.proposal,
                       njobs = 90,
                       tmp_path = "/gpfs/scratch/dw2625",
                       job_name = "BSIM_55",
                       sbatch_opt = list(time = "12:00:00",partition = "cpu_short", `mem-per-cpu` = "25G"),
                       export = c("bsim.fn","Sbeta","logDen","vonmises.lpr","metrop","pred.fn",
                                  "generate.data2", "bayes_single_model", "pred_sng_bayes"),
                       plan = "wait",
                       overwrite=TRUE)
  res <- Slurm_collect(sjob) # data is a list
  #res<- site_plasma_all[lapply(site_plasma_all, function(x) length(x))>1] #filter out the error message
  res <- rbindlist(res) # converting list to data.table
  
  date_stamp <- gsub("-", "", Sys.Date())
  dir.create(file.path("/gpfs/home/dw2625/r/BSIM/", date_stamp), showWarnings = FALSE)
  save(res, file = paste0("/gpfs/home/dw2625/r/BSIM/", date_stamp, "/scenarios",r,".rda"))
  results.aggregated2[[r]] <- res
}

date_stamp <- gsub("-", "", Sys.Date())
dir.create(file.path("/gpfs/home/dw2625/r/BSIM/", date_stamp), showWarnings = FALSE)
save(results.aggregated2, file = paste0("/gpfs/home/dw2625/r/BSIM/", date_stamp, "/all_6_scenarios_spline.rda"))


