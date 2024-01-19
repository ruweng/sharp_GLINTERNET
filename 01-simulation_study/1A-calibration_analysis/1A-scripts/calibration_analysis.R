#### Simulation Study 01 ----

#' Research Question: 
#' is the most stable glinternet model, as defined by the max stability score,
#' performing the best in terms of identifying true predictors [variable selection] as depicted by the F1 score,
#' across different parameter combinations [λ, π] ?

#' Procedure: 
#' 1. Simulate Regression with pairwise interaction effects following strong hierarchy
#' 2. For each parameter combination [λ, π] compute stability score and calculate performance
#' 3. Plot performance as depicted by the F1-score as a function of the stability score

### Set Working Directory and Load Library ----

rm(list = ls())
path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

source("../../../00-utils/glinternet_utils.R")
source("../../../00-utils/glinternet_viz.R")
library(igraph)
library(ggplot2)
library(tidyverse)
library(pheatmap)

### Parameters ----

n = 1000
pk = 10
nu_xy = 0.3
nu_int = 0.5
ev_xy = 0.5
beta_abs = c(0.1,1)
beta_sign = c(-1,1)
hierarchy = "strong"
family = "gaussian"

### Simulate Regression with Pairwise Interaction Effects ----
  set.seed(1) # set seed for reproducibility
  sim_int <- SimulateInteraction(n = n, # observations
                                 pk = pk, # main effect variables
                                 nu_xy = nu_xy, # expected probability of true main effect
                                 nu_int = nu_int, # expected probability of interaction effect
                                 ev_xy = ev_xy, # expected proportion of variance explained by true predictors 
                                 beta_abs = beta_abs, # vector defining the range of nonzero regression coefficients
                                 beta_sign = beta_sign, # vector of possible signs for regression coefficients
                                 family = family, # regression model
                                 hierarchy = hierarchy # interactions only between true main effects
  )
  
  saveRDS(sim_int, file = "../1A-data/sim.rds")
  
  pdf(file = paste0("../1A-results/1A-figures/sim_cor.pdf"), width = 8, height = 8)
  pheatmap(cor(sim_int$xdata), col = colgradient(100))
  dev.off()
  
  xdata_int <- model.matrix(~ 0 + .^2, data = as.data.frame(sim_int$xdata))

  pdf(file = paste0("../1A-results/1A-figures/sim_int_cor.pdf"), width = 8, height = 8)
  pheatmap(cor(xdata_int), col = colgradient(100))
  dev.off()
  
### LASSO ----
  
### Perform Variable Selection via the cross-validated LASSO ----
  set.seed(1) 
  cv <- glmnet::cv.glmnet(xdata_int,
                          sim_int$ydata, 
                          nfolds = 10, 
                          intercept = FALSE)
  
  saveRDS(cv, file="../1A-results/1A-models/lasso_cv.rds")
  
  pdf(file = paste0("../1A-results/1A-figures/lasso_cvCalibration.pdf"), width = 8, height = 6)
  plot.cv_calibration(cv, model = "lasso")
  dev.off()
  
  pdf(file = paste0("../1A-results/1A-figures/lasso_cvRegPaths_lambdahat.pdf"), width = 12, height = 6)
  plot.reg_paths(reg_beta = NULL, cv = cv, calibration = "lambdahat", sim_int = sim_int)
  dev.off()
  
  pdf(file = paste0("../1A-results/1A-figures/lasso_cvRegPaths_lambdahat_sd.pdf"), width = 12, height = 6)
  plot.reg_paths(reg_beta = NULL, cv = cv, calibration = "lambdahat_sd", sim_int = sim_int)
  dev.off()
  
### Perform Variable Selection via the Stability selected LASSO ----

  ## perform stability selected lasso
  stab_z <- VariableSelection(xdata_int,
                              sim_int$ydata,
                              implementation = PenalisedRegression,
                              n_cat = NULL) # estimate stability score using z statistic 
  
  saveRDS(stab_z, file = "../1A-results/1A-models/lasso_stab_zscore.rds")
  
  stab_nll <- VariableSelection(xdata_int,
                                sim_int$ydata,
                                implementation = PenalisedRegression,
                                pi_list = seq(0.5, 0.99, by = 0.01),
                                n_cat = 3) # estimate stability score using negative log likelihood
  
  saveRDS(stab_nll, file = "../1A-results/1A-models/lasso_stab_nll.rds")
  
  ## plot stability selection calibration
  pdf(file = paste0("../1A-results/1A-figures/lasso_stabCalibration_zscore.pdf"), width = 8, height = 8)
  par(mar=c(8,8,8,8))
  CalibrationPlot(stab_z, col = c("#ffedc0", "#ff0000", "darkred"))
  dev.off()
  
  pdf(file = paste0("../1A-results/1A-figures/lasso_stabCalibration_nll.pdf"), width = 8, height = 8)
  par(mar=c(8,8,8,8))
  CalibrationPlot(stab_nll, col = c("#ffedc0", "#ff0000", "darkred"))
  dev.off()
  
  ## plot variable selection as selection proportion at max stability
  pdf(file = paste0("../1A-results/1A-figures/lasso_selpropPerformance_zscore.pdf"), width = 12, height = 4)
  par(mar=c(6,8,4,8))
  plot.selprop_performance(stab_z, theta = sim_int$theta)
  dev.off()
  
  pdf(file = paste0("../1A-results/1A-figures/lasso_selpropPerformance_nll.pdf"), width = 12, height = 4)
  par(mar=c(6,8,4,8))
  plot.selprop_performance(stab_nll, theta = sim_int$theta)
  dev.off()
  
  ## plot StabilityPaths - the paths of selection proportions (stability) across different lambdas 
  pdf(file = paste0("../1A-results/1A-figures/lasso_stabPaths_zscore.pdf"), width = 12, height = 6)
  plot.stab_paths(stab_z, sim_int)
  dev.off()
  
  pdf(file = paste0("../1A-results/1A-figures/lasso_stabPaths_nll.pdf"), width = 12, height = 6)
  plot.stab_paths(stab_nll, sim_int)
  dev.off()
  
  ## For each parameter combination [λ, π] compute Stability Score and calculate Model Performance [F1] ----
  stab_z_grid_performance <- grid_performance(stab_z,sim_int)
  saveRDS(stab_z_grid_performance, file = "../1A-results/1A-performance/lasso_stabGridPerformance_zscore.rds")
  
  pdf(file = paste0("../1A-results/1A-figures/lasso_stabGridPerformance_zscore.pdf"), width = 14, height = 8)
  plot.grid_performance(stab_z_grid_performance)
  dev.off()
  
  stab_nll_grid_performance <- grid_performance(stab_nll,sim_int)
  saveRDS(stab_nll_grid_performance, file = "../1A-results/1A-performance/lasso_stabGridPerformance_nll.rds")
  
  pdf(file = paste0("../1A-results/1A-figures/lasso_stabGridPerformance_nll.pdf"), width = 14, height = 8)
  plot.grid_performance(stab_nll_grid_performance)
  dev.off()
  
### GLINTERNET ----
  
### Perform Variable Selection via the cross-validated GLINTERNET ----
  
  ## calculate number of levels per variable
  numLevels <- apply(sim_int$xdata,2, function(x) {
    if (is.factor(x))
      length(levels(x))
    else {1}})
  
  ## perform cross-validated glinternet
  set.seed(1) 
  cv_int <- glinternet::glinternet.cv(X = sim_int$xdata,
                                      Y = sim_int$ydata,
                                      numLevels = numLevels,
                                      nFolds = 10)
  
  saveRDS(cv_int, file = "../1A-results/1A-models/glinternet_cv.rds")
  
  ## plot cross-validation calibration
  
  pdf(file = paste0("../1A-results/1A-figures/glinternet_cvCalibration.pdf"), width = 8, height = 6)
  plot.cv_calibration(cv_int, model = "glinternet")
  dev.off()
  
  ## plot regularization paths - betas as a function of lambda
  
  reg_int <- glinternet::glinternet(X = sim_int$xdata,
                                    Y = sim_int$ydata,
                                    numLevels = numLevels,
                                    lambda = cv_int$lambda)
            
  reg_beta <- glinternet.beta(reg_int, xdata = sim_int$xdata)
  
  pdf(file = paste0("../1A-results/1A-figures/glinternet_cvRegPaths_lambdahat.pdf"), width = 12, height = 6)
  plot.reg_paths(reg_beta = reg_beta, cv = cv_int, calibration = "lambdahat", sim_int = sim_int)
  dev.off()
  
  pdf(file = paste0("../1A-results/1A-figures/glinternet_cvRegPaths_lambdahat_sd.pdf"), width = 12, height = 6)
  plot.reg_paths(reg_beta = reg_beta, cv = cv_int, calibration = "lambdahat_sd", sim_int = sim_int)
  dev.off()

### Perform Variable Selection via the Stability selected GLINTERNET ----

  ## calculate a grid of penalty parameters using implementation glinternet
  Lambda_int <- LambdaGridRegression_glinternet(sim_int$xdata, sim_int$ydata)
  
  ## perform stability selected glinternet
  stab_int_z <- VariableSelection(sim_int$xdata,
                                sim_int$ydata,
                                implementation = PenalisedInteraction,
                                Lambda = Lambda_int,
                                n_cat = NULL) # estimate stability score using z statistic 
  
  saveRDS(stab_int_z, file = "../1A-results/1A-models/glinternet_stab_zscore.rds")
  
  stab_int_nll <- VariableSelection(sim_int$xdata,
                                  sim_int$ydata,
                                  implementation = PenalisedInteraction,
                                  Lambda = Lambda_int,
                                  pi_list = seq(0.5, 0.99, by = 0.01),
                                  n_cat = 3) # estimate stability score using negative log likelihood

  saveRDS(stab_int_nll, file = "../1A-results/1A-models/glinternet_stab_nll.rds")
  
  ## plot stability selection calibration
  pdf(file = paste0("../1A-results/1A-figures/glinternet_stabCalibration_zscore.pdf"), width = 8, height = 8)
  par(mar=c(8,8,8,8))
  CalibrationPlot(stab_int_z, col = c("#ffedc0", "#ff0000", "darkred"))
  dev.off()
  
  pdf(file = paste0("../1A-results/1A-figures/glinternet_stabCalibration_nll.pdf"), width = 8, height = 8)
  par(mar=c(8,8,8,8))
  CalibrationPlot(stab_int_nll, col = c("#ffedc0", "#ff0000", "darkred"))
  dev.off()
  
  ## plot variable selection as selection proportion at max stability
  pdf(file = paste0("../1A-results/1A-figures/glinternet_selpropPerformance_zscore.pdf"), width = 12, height = 4)
  par(mar=c(6,8,4,8))
  plot.selprop_performance(stab_int_z, theta = sim_int$theta)
  dev.off()
  
  pdf(file = paste0("../1A-results/1A-figures/glinternet_selpropPerformance_nll.pdf"), width = 12, height = 4)
  par(mar=c(6,8,4,8))
  plot.selprop_performance(stab_int_nll, theta = sim_int$theta)
  dev.off()
  
  ## plot StabilityPaths - the paths of selection proportions (stability) across different lambdas 
  pdf(file = paste0("../1A-results/1A-figures/glinternet_stabPaths_zscore.pdf"), width = 12, height = 6)
  plot.stab_paths(stab_int_z, sim_int)
  dev.off()
  
  pdf(file = paste0("../1A-results/1A-figures/glinternet_stabPaths_nll.pdf"), width = 12, height = 6)
  plot.stab_paths(stab_int_nll, sim_int)
  dev.off()

  ## For each parameter combination [λ, π] compute Stability Score and calculate Model Performance [F1] ----
  stab_z_grid_performance <- grid_performance(stab_int_z,sim_int)
  saveRDS(stab_z_grid_performance, file = "../1A-results/1A-performance/glinternet_stabGridPerformance_zscore.rds")
  
  pdf(file = paste0("../1A-results/1A-figures/glinternet_stabGridPerformance_zscore.pdf"), width = 14, height = 8)
  plot.grid_performance(stab_z_grid_performance)
  dev.off()
  
  stab_nll_grid_performance <- grid_performance(stab_int_nll,sim_int)
  saveRDS(stab_nll_grid_performance, file = "../1A-results/1A-performance/glinternet_stabGridPerformance_nll.rds")
  
  pdf(file = paste0("../1A-results/1A-figures/glinternet_stabGridPerformance_nll.pdf"), width = 14, height = 8)
  plot.grid_performance(stab_nll_grid_performance)
  dev.off()
  
  
  
  
  
 

  