## Bayes single-index models for binomial response 
## main function 
bsim.fn <- function(y,  # target data series (a vector of size n) 
                    A,  # treatment data series (a vector of size n) 
                    X,  # the matrix of predictors to be included in the treatment interaction, ie. an n x p matrix (n data points and p predictors)
                    Xm = NULL, # the matrix of predictors to be included in the main effect, ie. an n x q matrix.
                    family = "binomial", 
                    beta.prior = NULL,     # prior beta 
                    lambda.prior = 100,    # typically ranges from 100-700; 100 means weak belief on beta.prior.
                    lambda.proposal= 100,  # scale parameter used in the proposal distribution 
                    m.prior = NULL, 
                    Q= NULL,   # prior variance for the main effect coefficient 
                    beta.ini = NULL,  # initial value for the frequentist estimate
                    n.burning=10,   # number of burn-in
                    n.samples=200,  # number of actual samples to be collected
                    k = 6)      # dimension of spline basis for the function g. 
#rho = NULL,   # tuning parameter that determines the smoothness of g.  
#bd.eps = NULL, 
#sample_seed=1111)
{
  
  X <- as.matrix(X)
  Xm <- as.matrix(Xm)
  y <- as.vector(y)
  A <- as.vector(A)
  p <- ncol(X) #number of covariates for the interaction term
  n <- nrow(X)
  if(is.null(Xm)) Xm <- X
  q <- ncol(Xm) #number of covariates for the main effects term
  
  A <- as.numeric(as.factor(A))
  A.unique <- unique(sort(A))
  L <- length(A.unique)
  A.mean <- mean(A)
  Ac <<- A - A.mean
  AXm <- cbind(1, A=Ac, Xm) 
  dat <- data.frame(y, A, X=X, Xm=Xm)
  
  # initialize the parameters using a non-Bayesian single-index model (Park et al 2021 Biometrics)
  b.ini <- simml::simml(y, A, X, Xm, family=family, k=k, beta.ini=beta.ini) 
  beta.ini = betax  <- b.ini$beta.coef
  x.beta <- X %*% betax  
  
  m.fit <- gam(y~ AXm-1, family=family)
  offset.mx <- as.matrix(predict(m.fit))
  
  D <- splines::ns(x.beta, df = k)  
  Pen <- diff(diag(k), differences = 2)   #Construct the penalty matrix
  Pen2 <- t(Pen)%*% Pen 
  
  knots <- attr(D, "knots")   # specification of the knot location for the component function g.
  Boundary.knots <- attr(D, "Boundary.knots")
  
  A.indicator = matrix(0, n, L)
  for (a in 1:L) A.indicator[, a] = as.numeric(A == A.unique[a])
  pr.A = summary(as.factor(A))/n
  C = NULL
  for (a in 1:L) C = cbind(C, diag(rep(pr.A[a], k)))
  qrc = qr(t(C))
  null.space.C = qr.Q(qrc, complete = TRUE)[, (nrow(C) + 1):ncol(C)] 
  
  # prior beta; if not given, we may use a frequentist estimate (and apply a diffuse prior variance)  
  if(is.null(beta.prior))  beta.prior <- betax 
  if(is.null(m.prior)) m.prior <- coef(m.fit)  
  if(is.null(Q)) Q <- diag(q+2)*10  # diffuse prior
  
  #Pen2 <- b$Pen2  # k x k penalty matrix associated with the spline fit 
  #if(is.null(rho))  rho <- b$sp   # roughness penalty (rho) associated with the spline fit  
  
  # setting up variables
  beta.sample <- matrix(0, p, n.samples) 
  gamma.sample <- matrix(0, k, n.samples)
  m.sample <- matrix(0, q+2, n.samples)
  
  accept <- rep(0, n.samples+n.burning)
  acceptsample <- rep(0, n.samples)
  reps <- n.burning+n.samples
  
  for(i in 1:reps){
    print(i)
    ## choose a new seed for each sample
    #set.seed(sample_seed+i)
    
    # ONE SWEEP OF THE METROPOLIS-HASTINGS SAMPLER
    
    # sample a new beta from von Misses-Fisher
    betap <- as.vector(rmovMF(1, lambda.proposal*betax))  # lambda.proposal will determine the acceptance probability for betap.
    cor(betap, betax)
    if(betap[1] < 0) betap <- -betap
    
    # compute the acceptance probability 
    lpr.obj <- vonmises.lpr(dat, offset.mx, 
                            betap, betax, 
                            lambda.prior, beta.prior, 
                            k, family, knots, Boundary.knots, Pen2, 
                            null.space.C, A.indicator)
    
    metrop.obj <- metrop(lpr.obj$lpr, betap, betax)  
    metrop.obj$accept
    accept[i]  <- metrop.obj$accept
    accept[i] 
    
    # update beta 
    betax  <- metrop.obj$keepv  
    if(accept[i]){
      Sbeta.obj <- lpr.obj$Sbetap
    }else{
      Sbeta.obj <- lpr.obj$Sbetax
    }    
    
    # update gamma 
    mu.gamma <- Sbeta.obj$gamma.hat
    Sigma.gamma <- Sbeta.obj$Sigma0/2
    gammax <- MASS::mvrnorm(1, mu.gamma, Sigma.gamma)
    
    offset.g <- Sbeta.obj$Dc %*% gammax
    
    # update main.effect 
    # need to compute W.check and z.check to get conditional posterior for m.  
    b <- gam(y~ offset(offset.g)+ AXm-1,  family=family)
    eta.hat <- b$offset + predict(b)
    mu.hat  <- b$family$linkinv(eta.hat)
    z <- (y - mu.hat)/ b$family$mu.eta(eta.hat)  + eta.hat 
    w <- b$family$mu.eta(eta.hat)^2/b$family$variance(mu.hat)
    W <- diag(w)
    Sigma.m <- solve(solve(Q) + t(AXm) %*% W %*% AXm)
    mu.m <- Sigma.m %*%(solve(Q)%*%m.prior + t(AXm) %*% W %*% z)
    mx <- MASS::mvrnorm(1, mu.m , Sigma.m)
    offset.mx <- predict(b, type="lpmatrix") %*% mx 
    
    
    # We only keep what happens after the burning period
    if(i > n.burning){
      j = i-n.burning
      acceptsample[i-n.burning] <- accept[i]
      m.sample[,j] = mx
      beta.sample[,j] = betax
      gamma.sample[,j] = gammax 
    }
    
  } # End of MH sampler
  
  
  list(beta.sample=beta.sample, 
       gamma.sample=gamma.sample, 
       m.sample=m.sample, 
       acceptsample=acceptsample, 
       knots=knots,  
       Boundary.knots=Boundary.knots, k=k, Pen2=Pen2, #rho=rho, 
       beta.prior=beta.prior, beta.ini=beta.ini, b.ini=b.ini, family=family, 
       null.space.C =null.space.C, A.mean=A.mean, A.unique=A.unique, L=L, y=y, X=X, Xm=Xm, A=A)
}






# computation of the quadratic form S(beta) with parameters 
# the design matrix X, the vector beta, the vector y 
Sbeta <- function(dat, offset.mx,  
                  betax, k, family, 
                  knots, Boundary.knots, 
                  Pen2, null.space.C, A.indicator)
{
  
  X  <- dat[,grep("X.", names(dat), fixed = TRUE)]
  A  <- dat[,2]
  Xm <- dat[,grep("Xm", names(dat))]
  y  <- dat[,1]
  
  # compute the index. 
  x.beta <-  as.matrix(X) %*% betax
  
  tmp0 = splines::ns(x.beta, df=k, knots=knots, Boundary.knots=Boundary.knots)
  tmp1 = NULL
  for (a in 1:ncol(A.indicator)) tmp1 = cbind(tmp1, A.indicator[, a]*tmp0)
  D  = tmp1
  Dc = D %*% null.space.C 
  
  # now we need to compute W and z to get conditional posterior for beta. 
  b <- gam(y~ offset(offset.mx)+ Dc-1,  family=family, paraPen=list(Dc =list(Pen2)))
  gamma.hat <- coef(b) # this is an adjusted coefficient (point estimate) adjusted for the effect of offset. 
  eta.hat <- b$offset + predict(b)
  mu.hat  <- b$family$linkinv(eta.hat)
  z <- (y - mu.hat)/ b$family$mu.eta(eta.hat)  + eta.hat 
  w <- b$family$mu.eta(eta.hat)^2/b$family$variance(mu.hat)
  W <- diag(w)
  
  aux1 <- t(Dc)%*% W %*% Dc
  aux2 <- t(Dc)%*% W %*% z
  if(det(aux1) < 10^{-2}){
    print("computationally singular")
    aux1 <- 0.01*diag(k) + t(Dc)%*% W %*% Dc  
  }
  Sigma0 <- solve(aux1)
  Sigmarho <- solve(aux1 + b$sp*Pen2)
  S1 <- t(z)%*%W%*%z + t(aux2) %*% Sigmarho %*% aux1 %*% Sigmarho %*% aux2
  aux3 <- diag(k) + Sigmarho %*% aux1
  Lambda <- aux3 %*% Sigma0 %*% t(aux3)
  
  S2 <- t(aux2) %*% Lambda %*% aux2 /4
  #tmp <- 0.5*Sigma0 %*% (diag(dim(aux1)[1]) + aux1 %*%  Sigmarho) %*% aux2 
  
  list(S1=S1, S2=S2, Sigma0=Sigma0, Sigmarho=Sigmarho, 
       Dc=Dc, 
       betax=betax, 
       offset.mx=offset.mx, 
       gamma.hat=gamma.hat, 
       eta.hat = eta.hat, 
       aux1=aux1, aux2=aux2, aux3=aux3)  
}



# Den - Integrated Density with respect to sigma and beta
logDen <- function(dat, offset.mx, 
                   betax, k, family, 
                   knots, Boundary.knots, Pen2, 
                   null.space.C, A.indicator)
{
  Sbeta.obj <- Sbeta(dat, offset.mx, 
                     betax, k, family, 
                     knots, Boundary.knots, Pen2, null.space.C, A.indicator) 
  logDen <- Sbeta.obj$S2  - 0.5*Sbeta.obj$S1
  
  list(logDen=logDen, Sbeta.obj=Sbeta.obj)
}




# VONMISESLPR - VON MISES Distribution - Log Probability Ratio
# returns the log of the D(betap) / D(betax) when both 
# betap and betax are distributed as Vonmises(lambda).
vonmises.lpr <- function(dat, offset.mx, 
                         betap, betax, 
                         lambda.prior, beta.prior, 
                         k, family, knots, Boundary.knots, Pen2, 
                         null.space.C, A.indicator)
{
  tmp1 <- logDen(dat, offset.mx, betap, k, family, 
                 knots, Boundary.knots, Pen2, null.space.C, A.indicator)
  tmp2 <- logDen(dat, offset.mx, betax, k, family, 
                 knots, Boundary.knots, Pen2, null.space.C, A.indicator)
  
  lpr <- tmp1$logDen - tmp2$logDen + lambda.prior* t(betap - betax) %*% beta.prior  
  
  list(lpr=lpr, Sbetap=tmp1$Sbeta.obj, Sbetax=tmp2$Sbeta.obj)
}


# METROP - perform a Metropolis-Hastings step
# LOGQ = log (unnormalized) density ratio
# If targ(x) is the target density (may be unnormalized)
# and gen(x) is the generating density (may be unnormalized),
# then LOGQ = log( targ(NV) gen(OV|NV) / targ(OV) gen(NV|OV) ) 
# In our case gen is symetric so
# LOGQ = log( targ(NV) / targ(OV)  ) 
# NV = new sample value
# OV = old sample value
# KEEPV = value of NV or OV that is kept
# ACCEPT = 1 if NV was kept, 0 if OV was kept
metrop <- function(logq, newval, oldval){
  if (logq>0){
    keepv = newval 
    accept = 1 
  }else{
    if(runif(1) < exp(logq)){
      keepv = newval
      accept = 1
    }else{
      keepv = oldval
      accept  = 0
    }
  }
  return(list(keepv=keepv, accept=accept))
}



# xnew: n by p matrix 
# this function generates a n by n.sample matrix, where each row records the distribution of each individual. 
pred.fn <- function(bsim.obj, xnew=NULL, xmnew=NULL, newA=NULL, thresh=0){
  
  if(is.null(xnew))  xnew  <- bsim.obj$X
  if(is.null(xmnew)) xmnew <- bsim.obj$Xm
  if(is.null(newA))  newA  <- bsim.obj$A
  
  xnew <- as.matrix(xnew)
  xmnew <- as.matrix(xmnew)
  newA <- as.numeric(as.factor(newA))  
  Axmnew <- cbind(1, newA-bsim.obj$A.mean, xmnew)
  
  null.space.C <- bsim.obj$null.space.C 
  n <- nrow(xnew)
  p <- ncol(xnew)
  q <- ncol(xmnew)
  
  beta.sample <- bsim.obj$beta.sample
  gamma.sample <- bsim.obj$gamma.sample
  m.sample <- bsim.obj$m.sample
  nsample <-  ncol(beta.sample)
  
  index.sample <- xnew %*% beta.sample   # n by nsample matrix
  eta.distr = contrast.distr = tbi.tmp <- matrix(rep(0, n*nsample), n, nsample)
  
  A.indicator = matrix(0, n, bsim.obj$L)
  for (a in 1:bsim.obj$L)  A.indicator[, a] = as.numeric(newA == bsim.obj$A.unique[a])
  
  for(i in 1:n){
    D.tmp1.i <- splines::ns(index.sample[i, ], df = bsim.obj$k, 
                            knots = bsim.obj$knots, Boundary.knots = bsim.obj$Boundary.knots) 
    D.tmp2.i = NULL
    for(a in 1:bsim.obj$L){
      D.tmp2.i = cbind(D.tmp2.i, D.tmp1.i * A.indicator[i, a])
    }
    gamma.tilde.sample <- null.space.C %*% gamma.sample
    
    eta.distr[i, ] <- as.matrix(t(Axmnew[i,])) %*% m.sample + rowSums(D.tmp2.i * t(gamma.tilde.sample))  
    contrast.distr[i, ] <- m.sample[2, ] + D.tmp1.i %*% (gamma.tilde.sample[bsim.obj$k+1:bsim.obj$k]- gamma.tilde.sample[1:bsim.obj$k])
    tbi.tmp[i, ] <- (contrast.distr[i,] < thresh)
  }
  
  #if(type=="response"){
  #  res <- bsim.obj$family$linkinv(eta.distr)
  #}
  tbi <- apply(tbi.tmp, 1, mean)
  
  list(eta.distr=eta.distr, contrast.distr =contrast.distr, tbi=tbi, index.sample=index.sample)
}
