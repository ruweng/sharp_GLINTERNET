#### Simulation Study 01 C Asymptotic Analysis ----

#' Research Questions:
#' Does the stability selected GLINTERNET perform better (on average across K iterations)
#' in identifying true predictors [Variable Selection] as depicted by the F1 score
#' than the cross-validated LASSO, the stability selected LASSO, and cross-validated GLINTERNET [standard case],
#' taking into account;
#' 1. differing dimensions [n x p] 
#' 2. differing proportion of effect sizes between main and interaction effects [beta_abs_int]
#' 
#' Procedure:
#' The analysis is run via an array job,
#' where each iteration constitutes its own job,
#' with parallelised models run sequentially.

### Load Libraries and Define Parameters ----

## Libraries

source("../../../../00-utils/glinternet_utils.R")
source("../../../../00-utils/glinternet_viz.R")

## PBS Arguments 

args = commandArgs(trailingOnly=TRUE)
n_cores = as.numeric(args[1])
node = as.numeric(args[2]) + 1

## Parameters

n <- 1000 # observations
pk <- c(50, 100) # number of variables
beta_abs_int <- seq(1, 0.1, by = -0.1)# expected proportion of explained variance
k <- 1:100 # iterations

## Parameter combinations 

params_grid <- expand.grid(k, beta_abs_int, pk, n)
colnames(params_grid) <- c("k","beta_abs_int", "pk", "n")

#### Simulations ----

simulate_data <- function(params){
  
  # Define parameters
  k_tmp <- params$k
  n_tmp <- params$n
  pk_tmp <- params$pk
  beta_abs_int_tmp <- params$beta_abs_int
  
  # Set seed based on iteration k
  set.seed(k_tmp)
  
  # Simulate data
  sim_int <- SimulateInteraction(n = n_tmp, # number of observations
                                 pk = pk_tmp, # number of variables
                                 ev_xy = 0.25, # expected proportion of variance explained by true predictors 
                                 nu_xy = 0.1, # expected probability of true main effect
                                 nu_int = 0.2, # expected probability of interaction effect
                                 hierarchy = "strong", # interactions only between true main effects
                                 family = "gaussian", # regression model
                                 beta_abs = 1, # vector defining the range of nonzero main effect regression coefficients
                                 beta_abs_int = beta_abs_int_tmp, # vector defining range of non-zero interaction effect regression coefficients
                                 beta_sign = 1) # vector of possible signs for regression coefficients
  return(sim_int)
}

#### GLINTERNET ----

## GLINTERNET Stability Selection

perform_stab_int <- function(sim_int, n_cat, n_cores, pi_list){
  
  Lambda_int <- LambdaGridRegression_glinternet(sim_int$xdata, sim_int$ydata)
  
  stab_int <- VariableSelection(sim_int$xdata,
                                sim_int$ydata,
                                implementation = PenalisedInteraction,
                                Lambda = Lambda_int,
                                tau = 0.5,
                                K = 100,
                                n_cat = n_cat, 
                                n_cores = n_cores, 
                                pi_list = pi_list)
  
  return(stab_int)
  
}

## GLINTERNET Cross-validation

perform_cv_int <- function(sim_int){
  
  numLevels <- apply(sim_int$xdata, 2, function(x){
    if (is.factor(x))
      length(levels(x))
    else 1})  
  
  cv_int <- glinternet::glinternet.cv(sim_int$xdata,
                                      sim_int$ydata,
                                      numLevels, 
                                      nFolds = 10)
  return(cv_int)
}

## GLINTERNET Retrieve theta

# Retrieve cv_int_theta at lambdahat
get_cv_int_theta_lambdahat <- function(cv_int, sim_int){
  
  numLevels <- apply(sim_int$xdata, 2, function(x){
    if (is.factor(x))
      length(levels(x))
    else 1})  
  
  lambdahat_tmp <- cv_int$lambdaHat
  
  glinternet_lambdahat <- glinternet::glinternet(sim_int$xdata, sim_int$ydata, lambda = lambdahat_tmp, numLevels = numLevels)
  
  theta_lambdahat <- glinternet.theta(mymodel = glinternet_lambdahat, xdata = sim_int$xdata)
  
  return(theta_lambdahat)
}

# Retrieve cv_int_theta at lambdahat_sd
get_cv_int_theta_lambdahat_sd <- function(cv_int, sim_int){
  
  numLevels <- apply(sim_int$xdata, 2, function(x){
    if (is.factor(x))
      length(levels(x))
    else 1})  
  
  lambdahat_sd_tmp <- cv_int$lambdaHat1Std
  
  glinternet_lambdahat_sd<- glinternet::glinternet(sim_int$xdata, sim_int$ydata, lambda = lambdahat_sd_tmp, numLevels = numLevels)
  
  theta_lambdahat_sd <- glinternet.theta(mymodel = glinternet_lambdahat_sd, xdata = sim_int$xdata)
  
  return(theta_lambdahat_sd)
}

#### LASSO ----

## LASSO Stability Selection

perform_stab <- function(sim_int, n_cat, n_cores, pi_list){
  
  xdata_int <- model.matrix(~ 0 + .^2, data = as.data.frame(sim_int$xdata))
  
  stab <- VariableSelection(xdata = xdata_int,
                            ydata = sim_int$ydata,
                            implementation = PenalisedRegression,
                            tau = 0.5,
                            K = 100,
                            n_cat = n_cat, 
                            n_cores = n_cores, 
                            pi_list = pi_list)
  
  return(stab)
}

## LASSO cross-validation

perform_cv <- function(sim_int){
  
  xdata_int <- model.matrix(~ 0 + .^2, data = as.data.frame(sim_int$xdata))
  
  cv <- glmnet::cv.glmnet(xdata_int,
                          sim_int$ydata, 
                          nfolds = 10)
  
  return(cv)
}

get_cv_theta_lambdahat <- function(cv, sim_int){
  
  xdata_int <- model.matrix(~ 0 + .^2, data = as.data.frame(sim_int$xdata))
  
  lambdahat_tmp <- cv$lambda.min
  
  lasso_lambdahat <- glmnet::glmnet(xdata_int,
                                    sim_int$ydata,
                                    lambda = lambdahat_tmp)
  
  coef_lambdahat <- as.matrix(coef(lasso_lambdahat)[-1, ,drop = F])
  
  theta_lambdahat <- replace(coef_lambdahat, which(coef_lambdahat != 0), 1)
  
  return(theta_lambdahat)
}

get_cv_theta_lambdahat_se <- function(cv, sim_int){
  
  xdata_int <- model.matrix(~ 0 + .^2, data = as.data.frame(sim_int$xdata))
  
  lambdahat_se_tmp <- cv$lambda.1se
  
  lasso_lambdahat_se <- glmnet::glmnet(xdata_int,
                                       sim_int$ydata,
                                       lambda = lambdahat_se_tmp)
  
  coef_lambdahat_se <- as.matrix(coef(lasso_lambdahat_se)[-1, ,drop = F])
  
  theta_lambdahat_se <- replace(coef_lambdahat_se, which(coef_lambdahat_se != 0), 1)
  
  return(theta_lambdahat_se)
}

#### Performance ----

get_performance <- function(sim_int,
                            stab_theta_z,
                            stab_theta_nll,
                            cv_theta_lambdahat,
                            cv_theta_lambdahat_se,
                            stab_int_theta_nll,
                            stab_int_theta_z,
                            cv_int_theta_lambdahat,
                            cv_int_theta_lambdahat_sd){
  
  # Retrieve model performance
  
  # LASSO
  performance_stab_z <- SelectionPerformance(theta = stab_theta_z, theta_star = sim_int$theta)
  performance_stab_nll <- SelectionPerformance(theta = stab_theta_nll, theta_star = sim_int$theta)
  performance_cv_lambdahat <- SelectionPerformance(theta = cv_theta_lambdahat, theta_star = sim_int$theta)
  performance_cv_lambdahat_se <- SelectionPerformance(theta = cv_theta_lambdahat_se, theta_star = sim_int$theta)
  
  # GLINTERNET
  performance_stab_int_z <- SelectionPerformance(theta = stab_int_theta_z, theta_star = sim_int$theta)
  performance_stab_int_nll <- SelectionPerformance(theta = stab_int_theta_nll, theta_star = sim_int$theta)
  performance_cv_int_lambdahat <- SelectionPerformance(theta = cv_int_theta_lambdahat, theta_star = sim_int$theta)
  performance_cv_int_lambdahat_sd <- SelectionPerformance(theta = cv_int_theta_lambdahat_sd, theta_star = sim_int$theta)
  
  # Add model name
  
  # LASSO
  performance_stab_z$model <- "lasso_stab_z"
  performance_stab_nll$model <- "lasso_stab_nll"
  performance_cv_lambdahat$model <- "lasso_cv_lambdahat"
  performance_cv_lambdahat_se$model <- "lasso_cv_lambdahat_se"
  
  # GLINTERNET
  performance_stab_int_z$model <- "glinternet_stab_z"
  performance_stab_int_nll$model <- "glinternet_stab_nll"
  performance_cv_int_lambdahat$model <- "glinternet_cv_lambdahat"
  performance_cv_int_lambdahat_sd$model <- "glinternet_cv_lambdahat_sd"
  
  # Combine objects
  performance <- rbind(performance_stab_z,
                       performance_stab_nll,
                       performance_cv_lambdahat,
                       performance_cv_lambdahat_se,
                       performance_stab_int_z, 
                       performance_stab_int_nll,
                       performance_cv_int_lambdahat,
                       performance_cv_int_lambdahat_sd)
  
  return(performance)
}    

#### Analysis ----

iterate <- function(params_grid, n_cores, node){
  
  # Define parameters
  params <- params_grid[node,]
  k_tmp <- params[,1]
  beta_abs_int_tmp <- params[,2]
  pk_tmp <- params[,3]
  n_tmp <- params[,4]
  
  # Run simulations
  sim_int <- simulate_data(params = params)
  saveRDS(sim_int, file = paste0("../1C-iii-data/sim_",node,".rds"))
  
  # Perform models
  
  # LASSO
  stab_z <- perform_stab(sim_int = sim_int, n_cat = NULL, n_cores = n_cores, pi_list = seq(0.01, 0.99, by = 0.01))
  saveRDS(stab_z, file= paste0("../1C-iii-results/1C-iii-models/lasso_stab_zscore_",node,".rds"))
  
  stab_nll <- perform_stab(sim_int = sim_int, n_cat = 3, n_cores = n_cores, pi_list = seq(0.5, 0.99, by = 0.01))
  saveRDS(stab_nll, file= paste0("../1C-iii-results/1C-iii-models/lasso_stab_nll_",node,".rds"))
  
  cv <- perform_cv(sim_int = sim_int)
  saveRDS(cv, file= paste0("../1C-iii-results/1C-iii-models/lasso_cv_",node,".rds"))
  
  # GLINTERNET
  stab_int_z <- perform_stab_int(sim_int = sim_int, n_cat = NULL, n_cores = n_cores, pi_list = seq(0.01, 0.99, by = 0.01)) 
  saveRDS(stab_int_z, file= paste0("../1C-iii-results/1C-iii-models/glinternet_stab_zscore_",node,".rds"))
  
  stab_int_nll <- perform_stab_int(sim_int = sim_int, n_cat = 3, n_cores = n_cores, pi_list = seq(0.5, 0.99, by = 0.01))
  saveRDS(stab_int_nll, file= paste0("../1C-iii-results/1C-iii-models/glinternet_stab_nll_",node,".rds"))
  
  cv_int <- perform_cv_int(sim_int = sim_int)
  saveRDS(cv_int, file= paste0("../1C-iii-results/1C-iii-models/glinternet_cv_",node,".rds"))
  
  # Retrieve thetas 
  
  # LASSO
  stab_theta_z <- SelectedVariables(stab_z)
  stab_theta_nll <- SelectedVariables(stab_nll)
  cv_theta_lambdahat <- get_cv_theta_lambdahat(sim_int = sim_int, cv = cv)
  cv_theta_lambdahat_se <- get_cv_theta_lambdahat_se(sim_int = sim_int, cv = cv)
  
  # GLINTERNET
  stab_int_theta_z <- SelectedVariables(stab_int_z)
  stab_int_theta_nll <- SelectedVariables(stab_int_nll)
  cv_int_theta_lambdahat <- get_cv_int_theta_lambdahat(sim_int = sim_int, cv = cv)
  cv_int_theta_lambdahat_sd <- get_cv_int_theta_lambdahat_sd(sim_int = sim_int, cv = cv)
  
  # Calculate performance
  performance <- get_performance(sim_int = sim_int,
                                 stab_theta_z = stab_theta_z,
                                 stab_theta_nll = stab_theta_nll,
                                 cv_theta_lambdahat = cv_theta_lambdahat,
                                 cv_theta_lambdahat_se = cv_theta_lambdahat_se,
                                 stab_int_theta_z = stab_int_theta_z,
                                 stab_int_theta_nll = stab_int_theta_nll,
                                 cv_int_theta_lambdahat = cv_int_theta_lambdahat,
                                 cv_int_theta_lambdahat_sd = cv_int_theta_lambdahat_sd)
  
  # Add parameters
  performance$k <- k_tmp
  performance$beta_abs_int <- beta_abs_int_tmp
  performance$pk <- pk_tmp
  performance$n <- n_tmp
  rownames(performance) <- NULL
  
  saveRDS(performance, file= paste0("../1C-iii-results/1C-iii-performances/performance_",node,".rds"))
  
  return(performance)
}

#### Run Analysis

effect_size_analysis <- iterate(params_grid = params_grid,
                                n_cores = n_cores,
                                node = node)


