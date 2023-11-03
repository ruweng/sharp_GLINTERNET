#' Intermediate repository for utility functions to be implemented in sharp
#' 
#' 
#' 
#' 
#' 
#' 
 
PenalisedInteraction <- function(xdata,
                                  ydata,
                                  Lambda = NULL,
                                  family,
                                  ...){
  
  # Calculating number of levels per variable
  numLevels <- apply(xdata,2, function(x) {
    if (is.factor(x))
      length(levels(x))
    else {1}})
  
  # Retrieving names of categorical and continuous variables
  xdata <- as.data.frame(xdata)
  cat_names <- colnames(xdata[,sapply(xdata, is.factor)])
  level_names <- do.call(rbind, lapply(xdata[,cat_names], levels))
  cont_names <- colnames(xdata[,sapply(xdata, is.numeric)])
  
  # Ordering categorical and continuous variables
  xdata_cat <- xdata[,cat_names]
  xdata_cont <- xdata[,cont_names]
  xdata <- cbind(xdata_cat, xdata_cont)
  
  # Changing factor levels to fit the glinternet required format [numeric characters ranging from 0:(nlevels-1)]
  xdata[,cat_names] <- apply(xdata[,cat_names],2,  FUN = function(x){factor(x, labels  = c(0:(length(unique(x))-1)))})
  
  # Storing extra arguments
  extra_args <- list(...)
  
  # Extracting relevant extra arguments
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = glinternet::glinternet)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("X", "Y", "lambda", "numLevels", "family", "maxIter")]
  
  # Running model
  mymodel <- do.call(glinternet::glinternet, args = c(
    list(
      X = xdata,
      Y = ydata,
      lambda = Lambda,
      numLevels = numLevels,
      family = family,
    ),
    tmp_extra_args
  ))
  
  # Retrieving coefficients from model
  coef <- coef(mymodel)
  
  # Changing categorical variables back to factors
  xdata[,cat_names] <- lapply(xdata[,cat_names], as.factor)
  
  # Generating matrix for main effects
  selected_main <- matrix(0, ncol = sum(numLevels), nrow = length(Lambda)) 
  
  
  attr(selected_main, "dimnames")[[2]] <- attr(model.matrix(~ 0 + .,
                                                            data = as.data.frame(xdata), 
                                                            contrasts.arg =  lapply(xdata[,cat_names],
                                                                                    contrasts,
                                                                                    contrasts=FALSE)),"assign")
  
  # Retrieving indices of main effects
  index_main_cat <- mapply("[",mapply("[",coef)[1,])[1,] # index of categorical variables
  index_main_cont <-  lapply(mapply("[",mapply("[",coef)[1,])[2,], FUN = function(x){x +length(xdata[,cat_names])}) # index of continuous variables
  index_main <- mapply(FUN = function(x,y){append(x, y)}, x = index_main_cat, y = index_main_cont)
  
  # Recording main effects based on indices
  selected_main <- as.data.frame(do.call(rbind,lapply(index_main, FUN = function(x){ifelse(attr(selected_main, "dimnames")[[2]] %in% x, 1,0)})))
  
  # Retrieving betas for main effects
  beta_main_cat <- lapply(mapply("[",mapply("[",coef)[2,])[1,], unlist)
  beta_main_cont <- lapply(mapply("[",mapply("[",coef)[2,])[2,], unlist)
  beta_main <- mapply(FUN = function(x,y){append(x, y)}, x = beta_main_cat, y = beta_main_cont)
  
  # Recording betas for main effects based on indices
  beta_main <- do.call(rbind,apply(as.matrix(Lambda),1, FUN = function(x){
    i <- which(Lambda == x)
    beta_main <- replace(selected_main[i,], which(selected_main[i,] == 1), beta_main[[i]])
    return(beta_main)
  })
  )
  
  # Generating matrix for interaction effects
  matrix_interaction <- matrix(0, nrow = length(Lambda), ncol = ncol(model.matrix(~ 0 + .^2,
                                                                                  data = as.data.frame(xdata),
                                                                                  contrasts.arg = lapply(xdata[,cat_names],
                                                                                                         contrasts,
                                                                                                         contrasts=FALSE))))
  attr(matrix_interaction, "dimnames")[[2]] <- attr(model.matrix(~ 0 + .^2,
                                                                 data = as.data.frame(xdata), 
                                                                 contrasts.arg = lapply(xdata[,cat_names],
                                                                                        contrasts,
                                                                                        contrasts=F)), "assign")
  
  matrix_interaction <- matrix_interaction[, -which(attr(matrix_interaction, "dimnames")[[2]] %in% attr(model.matrix(~ 0 + .,
                                                                                                                     data = as.data.frame(xdata), 
                                                                                                                     contrasts.arg = lapply(xdata[,cat_names],
                                                                                                                                            contrasts,
                                                                                                                                            contrasts=FALSE)), "assign")), drop = FALSE] # keep only interactions
  
  # Creating a dictionary for interaction indices
  index_int_dict <- as.data.frame(cbind(combn(ncol(xdata), 2)[1,], combn(ncol(xdata), 2)[2,], unique(attr(matrix_interaction, "dimnames")[[2]])))
  colnames(index_int_dict) <- c("V1", "V2", "V1:V2")
  
  # Retrieving matrices of interaction indices
  index_int_catcat <- mapply("[",mapply("[",coef)[3,])[1,]
  index_int_contcont <- lapply(mapply("[",mapply("[",coef)[3,])[2,], FUN = function(x){x + length(xdata[,cat_names])})
  index_int_catcont <- lapply(mapply("[",mapply("[",coef)[3,])[3,], FUN = function(x){
    V1 <- x[,1]
    V2 <- x[,2] + length(xdata[,cat_names])
    out <- cbind(V1,V2)
    return(out)})
  
  index_int <- mapply(FUN = function(x,y,z){rbind(x,y,z)}, x = index_int_catcat, y = index_int_contcont, z = index_int_catcont)
  
  # Mapping model.matrix indices to dictionary of interaction indices
  index_int <- lapply(index_int, FUN = function(x){merge(x, index_int_dict, sort = F)})
  
  # Recording interaction effects based on mapped indices
  selected_int <- as.data.frame(do.call(rbind,lapply(index_int, FUN = function(x){
    idx <- x[,3]
    selected_int <- ifelse(attr(matrix_interaction, "dimnames")[[2]] %in% idx, 1,0)
    return(selected_int)
  })))
  
  attr(selected_int, "assign") <- attr(matrix_interaction, "dimnames")[[2]] 
  
  # Creating a dictionary for interaction betas
  beta_int_dict <- as.data.frame(attr(selected_int, "assign"))
  colnames(beta_int_dict) <- "V1:V2"
  beta_int_dict <- merge(index_int_dict, beta_int_dict)
  
  # Retrieving betas for interaction effects
  beta_int_catcat <- lapply(mapply("[",mapply("[",coef)[4,])[1,], unlist)
  beta_int_contcont <- lapply(mapply("[",mapply("[",coef)[4,])[2,], unlist)
  beta_int_catcont <- lapply(mapply("[",mapply("[",coef)[4,])[3,], unlist)
  
  beta_int_full <- mapply(FUN = function(x,y,z){as.vector(c(x,y,z))}, x = beta_int_catcat, y = beta_int_contcont, z = beta_int_catcont)
  
  # Assigning interaction indices to corresponding betas
  beta_index_int <- lapply(index_int, FUN = function(x){merge(x, beta_int_dict, sort = F)})
  
  beta_int_map <- mapply(FUN = function(idx,beta){
    beta_int_map <- cbind(idx,beta)},
    idx = beta_index_int,
    beta = beta_int_full)
  beta_int_map <- lapply(beta_int_map, FUN = function(x){x[order(as.numeric(x[,3])),]})
  
  # Recording betas for interaction effects based on mapped indices
  beta_int <- matrix_interaction
  beta_int <- t(mapply(FUN = function(x,y){
    i <- which(Lambda == x)
    idx <- y[,3]
    beta <- y$beta
    beta_int <- replace(beta_int[i,], which(colnames(beta_int) %in% idx), beta)
  }, x = as.matrix(Lambda), y = beta_int_map))
  
  # Combining matrices of main and interaction effects
  selected <- cbind(selected_main, selected_int)
  beta_full <- cbind(beta_main, beta_int)
  
  # Changing column and rownames
  colnames(selected) <- colnames(beta_full) <- colnames(model.matrix( ~ 0 + .^2,
                                                                      data = as.data.frame(xdata),
                                                                      contrasts.arg = lapply(X = lapply(xdata[,cat_names],
                                                                                                        FUN = function(x){
                                                                                                          constrasts <- contrasts(x, contrasts = F)}
                                                                      ),
                                                                      FUN = function(x){
                                                                        colnames(x) <- paste0("_", colnames(x))
                                                                        rownames(x) <- paste0("_", rownames(x))
                                                                        x
                                                                      })))
  
  rownames(selected) <- paste0("s",0:(length(Lambda)-1))
  rownames(beta_full) <- paste0("s",0:(length(Lambda)-1))
  
  # Returning output 
  selected <- as.matrix(selected)
  beta_full <- as.matrix(beta_full)
  
  return(list(selected = selected, beta_full = beta_full))
  
}



#' TODO add implementation argument to LambdaGridRegression to either call glinternet or glmnet



LambdaGridRegression_glinternet <- function(xdata, ydata, tau = 0.5, seed = 1,
                                            family = "gaussian",
                                            resampling = "subsampling",
                                            Lambda_cardinal = 100, check_input = TRUE,
                                            ...) {
  # Object preparation, error and warning messages
  Lambda <- NULL
  pi_list <- seq(0.6, 0.9, by = 0.01)
  K <- 100
  n_cat <- 3
  PFER_method <- "MB"
  PFER_thr <- Inf
  FDP_thr <- Inf
  verbose <- TRUE
  implementation <- PenalisedRegression
  # Checks are not re-run if coming from VariableSelection to avoid printing twice the same messages
  if (check_input) {
    # CheckInputRegression(
    #   xdata = xdata, ydata = ydata, Lambda = Lambda, pi_list = pi_list,
    #   K = K, tau = tau, seed = seed, n_cat = n_cat,
    #   family = family, implementation = implementation,
    #   resampling = resampling, PFER_method = PFER_method,
    #   PFER_thr = PFER_thr, FDP_thr = FDP_thr,
    #   Lambda_cardinal = Lambda_cardinal,
    #   verbose = verbose
    # )
    # Object preparation, error and warning messages
    CheckParamRegression(
      Lambda = Lambda, pi_list = pi_list,
      K = K, tau = tau, seed = seed, n_cat = n_cat,
      family = family, implementation = implementation,
      resampling = resampling, PFER_method = PFER_method,
      PFER_thr = PFER_thr, FDP_thr = FDP_thr,
      Lambda_cardinal = Lambda_cardinal,
      verbose = verbose
    )
    CheckDataRegression(
      xdata = xdata, ydata = ydata, family = family, verbose = verbose
    )
  }
  rm(n_cat)
  rm(Lambda)
  rm(pi_list)
  rm(K)
  
  numLevels <- apply(xdata,2, function(x) {
    if (is.factor(x))
      length(levels(x))
    else 1})
  
  # Taking one subsample/boostrap sample of the data
  withr::local_seed(1) # To keep to allow for reproducible parallelisation
  s <- Resample(data = ydata, family = family, tau = tau, resampling = resampling)
  
  # Applying function for variable selection to get upperbound of Lambda
  withr::local_seed(1) # To keep to allow for reproducible parallelisation
  mymodel <- glinternet::glinternet(X = xdata[s, ], Y = ydata[s, ], family = family, numLevels =  numLevels, nLambda = Lambda_cardinal)
  # mycv <- do.call(glmnet::glmnet, args = list(xdata = xdata[s, ], ydata = ydata[s, ], family = family, ...))
  
  # Creating a grid of lambda values from min and max
  Lambda <- cbind(LambdaSequence(lmax = max(mymodel$lambda), lmin = min(mymodel$lambda), cardinal = Lambda_cardinal))
  Lambda <- as.matrix(stats::na.exclude(Lambda))
  rownames(Lambda) <- paste0("s", seq(0, nrow(Lambda) - 1))
  
  return(Lambda)
} 