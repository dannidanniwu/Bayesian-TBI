# This function fit a simple Bayesian linear model for the outcome
bayes_spline_model <- function(X, p, y, A, mod, xnew=xnew, xmnew=xmnew, newA=newA, coefs = sgbayes.coefs,thresh = 0){
  #set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
  N = nrow(X)
  D =  ncol(X)
  J = length(unique(A))
  B = splines::bs(X, df = 5, intercept = TRUE)
  # matplot(dd$timeID, B, type = 'l', lty = 1, 
  #         xlab = "Time", ylab = "B-spline Basis", 
  #         main = "B-spline Basis Functions")
  
  K_bspline = k # Number of B-spline basis functions you desire (for example)
  # Initialize a list to store B-spline matrices for each covariate
  B_list <- vector("list", D)
  Bnew_list <- vector("list", D)
  knots <- list()
  Boundary.knots <- list()
  
  for (d in 1:D) {
    B_list[[d]] <- splines::ns(X[, d], df = K_bspline)
    knots[[d]] <- attr(B_list[[d]], "knots")   # specification of the knot location for the component function g.
    Boundary.knots[[d]] <- attr(B_list[[d]], "Boundary.knots")
  }
  
  Y <- y
  
  if(is.null(xnew))  xnew  <- bsim.obj$X
  if(is.null(xmnew)) xmnew <- bsim.obj$Xm
  if(is.null(newA))  newA  <- bsim.obj$A
  
  N_new = nrow(xnew)
  for (d in 1:D) {
    Bnew_list[[d]] <- splines::ns(xnew[, d], df = K_bspline, knots = knots[[d]], Boundary.knots = Boundary.knots[[d]])
  }
  
  studydata <- list(
    N = N,
    D = D,
    X = X,
    J = J, 
    K_bspline = K_bspline,
    B = B_list,
    Y = y, 
    A = A, 
    N_new = N_new,
    Xnew = xnew,
    Bnew = Bnew_list,
    Anew = newA)
  
  
  fit_spline <- mod$sample(
    data = studydata,
    refresh = 0,
    chains = 4L,
    parallel_chains = 4L,
    iter_warmup = 500,
    iter_sampling = 2500,
    show_messages = FALSE) 
  
  # Check model divergence and get posterior distributions of parameters 
  diagnostics_df <- as_draws_df(fit_spline$sampler_diagnostics())
  div_spline <- sum(diagnostics_df[, 'divergent__'])
  
  # Get posterior draws of all parameters
  draws_dt <- data.table(as_draws_df(fit_spline$draws()))
  
  eta_distr <- as.matrix(draws_dt%>%select(paste0("eta","[",1:N_new,"]"))) #N_sample*N_new
  contrast_distr <- as.matrix(draws_dt%>%select(paste0("contrast_distr","[",1:N_new,"]"))) #N_sample*N_new
  tbi <-  apply(contrast_distr, 2, function(column_mean) mean(column_mean < thresh))
  
  return(list(eta_distr=eta_distr, contrast_distr =contrast_distr, tbi=tbi, div=div_spline))
}

