library(simstudy)
library(data.table)
library(mgcv)
library(lme4)
library(splines)
library(brms)
library(cmdstanr)
library(slurmR)
library(dplyr)
# Compile the Stan model
mod <- cmdstan_model("./one_x_a_interaction_penalty.stan")
s_define <- function() {
  #x
  def <- defData(varname = "x1", formula = 0, variance = 1)
  def <- defData(def, varname = "x2", formula = 0, variance = 2)
  def <- defData(def, varname = "A", formula = "1;2", dist = "uniform")
  #y
  defOut <- defDataAdd(varname = "y", formula = "exp(-(x1 - 0.5 + A + x2)^2) + 5*A", variance = 0.25)
  #A=0: exp(-((x1+x2) - 0.5)^2) 
  #A=1: exp(-((2*x1+5*x2) - 0.5)^2) + 5
  #A=1: exp(-((2*(x1 + x2)) - 0.5)^2) + 5
  
  
  return(list(def = def, defOut = defOut)) 
}

s_generate <- function(iter, list_of_defs) {
  
  list2env(list_of_defs, envir = environment())
  set.seed(iter+123)
  #--- add data generation code ---#
  #24 sites in total
  ds <- genData(200, def)
  ds <- addColumns(defOut, ds)
  #summary(ds)
  ds[] #  generated_data is a data.table
}

s_model <- function(train_data, test_data, mod) {
  #set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
  #######Fitting the GAM model with default penalization 
  fitgam <- gam(y ~ A + s(x1, by=A,bs="bs",k=9+4) + s(x2, by=A,bs="bs",k=9+4) + s(x1, by=x2,bs="bs",k=9+4) -1, data = train_data, method="REML")
  res_fitgam <- c(summary(fitgam)$p.coeff["A"], summary(fitgam)$se["A"])
  range <-   res_fitgam[1] + c(-1,1) * 1.96 *   res_fitgam[2]
  
  # Define the knots based on the training data
  knots_train <- quantile(train_data$x1 * train_data$x2 * train_data$A, probs=seq(0, 1, length=11)[-c(1, 11)])
  
  B_train <- predict(bs(train_data$x1 * train_data$x2 * train_data$A, degree=3, knots=quantile(train_data$x1, probs=seq(0, 1, length=11)[-c(1, 11)])))
  colnames(B_train) <- paste0("Bspline_", 1:ncol(B_train))
  
  # Compute the B-spline basis for the test set using the same knots from the training data
  B_test_tr <- predict(bs(test_data$x1 * test_data$x2 * test_data$A, degree=3, knots=knots_train))
  colnames(B_test_tr) <- paste0("Bspline_", 1:ncol(B_test_tr))
  
  stan_data <- list(num_data = nrow(train_data),
                    num_basis = ncol(B_train),
                    B = t(B_train),
                    A = train_data$A,
                    y = train_data$y,
                    num_data_test = nrow(test_data),
                    B_test = t(B_test_tr),
                    A_test = test_data$A
  )
  
  # Fit the Bayesian model
  fit <- mod$sample(data = stan_data,
                    refresh = 0,
                    chains = 4L,
                    parallel_chains = 4L,
                    iter_warmup = 500,
                    iter_sampling = 2500,
                    show_messages = FALSE)
  
  diagnostics_df <- as_draws_df(fit$sampler_diagnostics())
  div <- sum(diagnostics_df[, 'divergent__'])
  bayes_gam = fit$summary(variables="beta_A",
                          posterior::default_summary_measures()[1:3],
                          quantiles = ~ quantile(., probs = c(0.025, 0.975)),
                          posterior::default_convergence_measures())
  covered_bayes =   (bayes_gam$`2.5%`< 5 & 5 < bayes_gam$`97.5%`)
  
  #Fit a frequentist linear model with the same basis as the Bayesian model, but no penalization
  # Incorporating the B-spline basis into the data
  ds_with_bspline <- cbind(train_data,  B_train)
  # Fitting the model using gam
  bspline_terms <- paste(colnames(B_train), collapse = " + ")
  formula_str <- paste("y ~ A +", bspline_terms, "- 1")
  fitgam2 <- gam(as.formula(formula_str), data = ds_with_bspline, method="REML")
  res_fitgam2 <- c(summary(fitgam2)$p.coeff["A"], summary(fitgam2)$se["A"])
  range2 <-   res_fitgam2[1] + c(-1,1) * 1.96 *   res_fitgam2[2]
  
  
  
  ##############Test RMSE on the test data###########gitgit######
  # Use the GAM model to predict on the test data to get the mean predictions.
  test_data$pred_gam <- predict(fitgam, newdata = test_data)
  
  
  #posterior_samples <- data.table(as_draws_df(fit$draws()))
  #a <- as.matrix(posterior_samples%>%select(paste0("a","[",1:length(knots_train)+3,"]")))
  #beta_A <- as.matrix(posterior_samples%>%select(paste0("beta_A")))
  
  posterior_predictions <-fit$summary(variables="y_pred_test",posterior::default_summary_measures()[1])
  test_data$pred_bayesian <- posterior_predictions$mean
  
  # Use the GAM nonpenalizaed version model to predict on the test data to get the mean predictions.
  ds_with_bspline_test <- cbind(test_data,B_test_tr)
  test_data$pred_gam_np <- predict(fitgam2, newdata = ds_with_bspline_test)
  
  # Calculate prediction error for all models
  test_data$error_gam <- (test_data$y - test_data$pred_gam)^2
  test_data$error_bayesian <- (test_data$y - test_data$pred_bayesian)^2
  test_data$error_gam_np <- (test_data$y - test_data$pred_gam_np)^2
  
  # Calculate RMSE for both models
  rmse_gam <- sqrt(mean(test_data$error_gam))
  rmse_bayesian <- sqrt(mean(test_data$error_bayesian))
  rmse_gam_np <- sqrt(mean(test_data$error_gam_np))
  
  model_results <- data.table(est_gam_freq= res_fitgam[1], se_gam_freq=res_fitgam[2], 
                              gam_lowci=range[1], gam_upci=range[2], bayes_gam,div, 
                              est_gam_freq_np = res_fitgam2[1], se_gam_freq_np=res_fitgam2[2], 
                              gam_np_lowci=range2[1], gam_np_upci=range2[2],
                              covered_gam_freq=(range[1] < 5 & 5 < range[2]),
                              covered_bayes,covered_gam_np_freq=(range2[1] < 5 & 5 < range2[2]),
                              rmse_gam, rmse_bayesian, rmse_gam_np) %>%
    mutate(across(-c(variable,covered_gam_freq, covered_bayes,covered_gam_np_freq), round, 3))
  
  setnames(model_results, c("est_gam_freq","se_gam_freq","lowci_freq", "upci_freq","variable","est_mean_bayes",
                            "est_med_bayes","est_sd_bayes",
                            "lowci_bayes",
                            "upci_bayes","rhat","ess_bulk","ess_tail",
                            "div","est_gam_np_freq","se_gam_np_freq","lowci_freq_np", "upci_freq_np",
                            "covered_freq","covered_bayes","covered_gam_np_freq",
                            "rmse_gam","rmse_bayesian","rmse_gam_np"))
  
  return(model_results)
}

s_single_rep <- function(iter,list_of_defs, mod) {
  
  train_data <- s_generate(iter,list_of_defs)
  test_data <- s_generate(iter+999,list_of_defs)
  model_results <- s_model(train_data, test_data, mod)
  
  return(model_results)
}

s_replicate <- function(iter, mod) {
  list_of_defs = s_define()
  model_results = s_single_rep(iter,list_of_defs, mod)
  return(data.table(iter=iter, model_results))
}


job <- lapply(1:10,function(i) s_replicate(iter=i, mod=mod))
res <- rbindlist(job)

res
sapply(res,mean)
