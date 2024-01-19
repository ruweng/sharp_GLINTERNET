#### GLINTERNET UTENSILS ####
# TODO: add plots for stability glinternet, i.e. IncrementalPlot (already possible!)
# TODO: add beta_int_abs?
# TODO: update LambdaGridRegression to be able to deal with glinternet
# TODO: clean up SampleInteraction and SimulateInteraction and make sure functionalities work
# TODO: implement into sharp and fake

### Simulate interactions ----

SampleInteractions <- function(xdata, theta, q = NULL,
                               nu_int = 0.1,
                               beta_abs = c(0.1, 1), beta_sign = c(-1, 1), continuous = TRUE, 
                               hierarchy = "strong"){
  
  # Definition of number of outcome variables
  if (is.null(q)) {
    q <- length(pk)
  }
  if (length(nu_int) != q) {
    nu_int <- rep(nu_int[1], q)
  }
  
  # Simulation of interactions
  n_theta <- vector(length = q) # number of variables per outcome
  n_pred <- vector(length = q) # number of true predictors per outcome
  n_noise <- vector(length = q) # number of noise variables per outcome
  
  pk <- vector(length = q, mode = "list") # block structure of interaction adjacency matrix
  nu_mat <- vector(length = q, mode = "list") # matrix of probability of interaction per outcome
  names <- vector(length = q, mode = "list") # list of sorted variable names per outcome
  
  int_adj <- vector(length = q, mode = "list") # interaction adjacency matrix per outcome
  int_table <- vector(length = q, mode = "list") # interaction table per outcome
  
  for (k in 1:q){
    
    # calculate adjacency matrix structure
    n_theta[k] <-  length(theta[,k])
    n_pred[k] <-  length(theta[,k][theta[,k]==1]) 
    n_noise[k] <-  length(theta[,k][theta[,k]==0])
    
    pk[[k]] <- c(n_pred[k], n_noise[k])
    pk[[k]] <- pk[[k]][pk[[k]] != 0] # remove zeros from pk
    
    if(length(pk[[k]]) == 1){
      if(hierarchy == "strong"){
        int_adj[[k]] <- SimulateAdjacency(pk = pk[[k]], nu_within = nu_int[k]) 
      } else if(hierarchy == "weak"){
        stop("There are no possible weak interactions. Please switch 'hierarchy 'to 'strong'")
      } else if(hierarchy == "anti"){
        stop("There are no possible non-hierarchical interactions. Please switch 'hierarchy' to 'strong'")
      } else("invalid input for argument 'hierarchy'. Must either be 'strong', 'weak' or 'anti'")
    } else{
      if(hierarchy == "strong"){
        nu_mat[[k]] <- matrix(c(nu_int[k], 0, 0, 0), nrow = 2)
      } else if(hierarchy == "weak"){
        nu_mat[[k]] <- matrix(c(0, nu_int[k], nu_int[k], 0), nrow = 2)
      } else if(hierarchy == "anti"){
        nu_mat[[k]] <- matrix(c(0, 0, 0, nu_int[k]), nrow = 2)
      } else{
        stop("invalid input for argument 'hierarchy'. Must either be 'strong', 'weak' or 'anti'")
      }
      # Simulate binary status for interactions with nu_int expected probability
      int_adj[[k]] <- SimulateAdjacency(pk = pk[[k]], nu_mat = nu_mat[[k]]) 
    }
    
    # add variable names 
    names[[k]] <- c(attr(theta[,k][theta[,k]==1], "names"), attr(theta[,k][theta[,k]==0], "names"))
    attr(int_adj[[k]], "dimnames") <- list(names[[k]], names[[k]])
    
    
    # convert into list of interactions per outcome
    int_adj[[k]][upper.tri(int_adj[[k]], diag = TRUE)] <- NA
    int_table[[k]] <- na.omit(data.frame(as.table(int_adj[[k]])))
    attr(int_table[[k]], "names") <- c("Var1", "Var2", "theta")
    
    
    # Sampling interaction coefficients if interactions are present
    int_table[[k]]["beta"] <- int_table[[k]]["theta"]
    if (continuous){
      int_table[[k]]["beta"] <- int_table[[k]]["beta"] * matrix(stats::runif(n = nrow(int_table[[k]]["beta"]), min = min(beta_abs), max = max(beta_abs)),
                                                                nrow = nrow(int_table[[k]]["beta"])
      )
    } else {
      int_table[[k]]["beta"] <- int_table[[k]]["beta"] * matrix(base::sample(beta_abs, size = nrow(int_table[[k]]["beta"]), replace = TRUE),
                                                                nrow = nrow(int_table[[k]]["beta"])
      )
    }
    int_table[[k]]["beta"] <- int_table[[k]]["beta"] * matrix(base::sample(beta_sign, size= nrow(int_table[[k]]["beta"]), replace = TRUE), 
                                                              nrow = nrow(int_table[[k]]["beta"])
    )
    # ordering interactions
    int_table[[k]]$Var1 <- ordered(int_table[[k]]$Var1, colnames(xdata))
    int_table[[k]]$Var2 <- ordered(int_table[[k]]$Var2, colnames(xdata))
    
    int_table[[k]] <- transform(int_table[[k]], Var1 = ifelse(Var1 > Var2, Var2, Var1), Var2 = ifelse(Var2 < Var1, Var1, Var2))
    int_table[[k]]$Var1 <- ordered(int_table[[k]]$Var1, 1:ncol(xdata))
    int_table[[k]]$Var2 <- ordered(int_table[[k]]$Var2, 1:ncol(xdata))
    
    levels(int_table[[k]]$Var1) <- colnames(xdata)
    levels(int_table[[k]]$Var2) <- colnames(xdata)
    
  }
  
  return(int_table) 
  
}

 SimulateInteraction <- function(n = 100, pk = 10, xdata = NULL,
                                family = "gaussian", q = 1,
                                theta = NULL, nu_xy = 0.2,
                                nu_int = 1, hierarchy = "strong",
                                beta_abs = c(0.1, 1), beta_sign = c(-1, 1), continuous = TRUE,
                                ev_xy = 0.7) {
  
  # TODO in future versions: introduce more families ("multinomial" and "cox")
  # Checking that either n and pk or xdata are provided
  if (is.null(xdata)) {
    if (is.null(pk) | is.null(n)) {
      stop("Argument 'xdata' must be provided if 'pk' and 'n' are not provided.")
    }
  }
  
  # Checking other inputs
  if (length(ev_xy) != q) {
    ev_xy <- rep(ev_xy[1], q)
  }
  if (length(nu_xy) != q) {
    nu_xy <- rep(nu_xy[1], q)
  }
  if (length(nu_int) !=q){
    nu_int <- rep(nu_int[1], q)
  }
  
  # Creating objects not provided as input
  if (!is.null(xdata)) {
    n <- nrow(xdata)
    p <- ncol(xdata)
  } else {
    p <- sum(pk)
    xsimul <- SimulateGraphical(
      n = n, pk = pk, theta = NULL,
      implementation = HugeAdjacency, topology = "random",
      nu_within = 0, nu_between = 0, nu_mat = NULL,
      v_within = 0, v_between = 0,
      v_sign = c(-1, 1), continuous = TRUE,
      pd_strategy = "diagonally_dominant", ev_xx = NULL, scale_ev = TRUE,
      u_list = c(1e-10, 1), tol = .Machine$double.eps^0.25,
      scale = TRUE, output_matrices = FALSE
    )
    xdata <- xsimul$data
  }
  
  # Checking theta if provided
  if (!is.null(theta)) {
    if (q == 1) {
      if (is.vector(theta)) {
        theta <- cbind(theta)
      }
    }
    if (ncol(theta) != q) {
      stop("Arguments 'theta' and 'q' are not compatible. Please provide a matrix 'theta' with 'q' columns.")
    }
    if (nrow(theta) != p) {
      stop("Please provide a matrix 'theta' with as many columns as predictors.")
    }
    theta <- ifelse(theta != 0, yes = 1, no = 0)
  }
  
  # Sampling true predictors
  if (is.null(theta)) {
    theta <- SamplePredictors(pk = p, q = q, nu = nu_xy, orthogonal = FALSE)
  }
  
  
  # Sampling regression coefficients
  beta <- theta
  if (continuous) {
    beta <- beta * matrix(stats::runif(n = nrow(beta) * ncol(beta), min = min(beta_abs), max = max(beta_abs)),
                          nrow = nrow(beta), ncol = ncol(beta)
    )
  } else {
    beta <- beta * matrix(base::sample(beta_abs, size = nrow(beta) * ncol(beta), replace = TRUE),
                          nrow = nrow(beta), ncol = ncol(beta)
    )
  }
  beta <- beta * matrix(base::sample(beta_sign, size = nrow(beta) * ncol(beta), replace = TRUE),
                        nrow = nrow(beta), ncol = ncol(beta)
  )
  
  # Define rownames of beta and theta
  rownames(beta) <- rownames(theta) <- colnames(xdata) 
  
  # Sampling true interactions (add option to implement own list of interactions and checks)
  
  # Sampling random network of pairwise interactions and interaction coefficients
  int <- SampleInteractions(xdata = xdata, theta = theta, q = q, nu_int = nu_int, hierarchy = hierarchy,
                            beta_abs = beta_abs, beta_sign = beta_sign, continuous = continuous)
  
  # Sampling outcome data
  ydata <- matrix(NA, ncol = q, nrow = nrow(xdata))
  if (family == "gaussian") {
    for (j in 1:q) {
      
      # Linear combination of main effects
      ypred_main <- xdata %*% beta[, j] 
      
      # Linear combination of interaction effects
      
        # matrix of xdata * interaction effects
        y_int <- matrix(0, ncol = nrow(int[[j]]), nrow = nrow(xdata)) 
        for (i in 1:nrow(int[[j]])){
          int_var1 <- as.character(int[[j]]$Var1[i])
          int_var2 <- as.character(int[[j]]$Var2[i])
          y_int[,i] <- int[[j]]$beta[i] * xdata[,int_var1] * xdata[,int_var2]
        }
        
      ypred_int <- rowSums(y_int)
      
      # Linear combination of effects
      ypred <- ypred_main + ypred_int
      
      # Calculating standard deviation that achieves desired proportion of explained variance
      sigma <- sqrt((1 - ev_xy[j]) / ev_xy[j] * stats::var(ypred)) # using var(ypred) is a shortcut (could use true variances but not always available)
      
      # Sampling from Normal distribution
      ydata[, j] <- stats::rnorm(n = n, mean = ypred, sd = sigma) # does sigma hold with interactions?
    }
  }
  if (family == "binomial") {
    for (j in 1:q) {
      # Linear combination of main effects
      main_crude_log_odds <- xdata %*% beta[, j]
      
      # Linear combination of interaction effects
      
      # matrix of xdata * interaction effects
      int_crude_log_odds_matrix <- matrix(NA, ncol = nrow(int[[j]]), nrow = nrow(xdata)) 
      for (i in 1:nrow(int[[j]])){
        int_var1 <- as.character(int[[j]]$Var1[i])
        int_var2 <- as.character(int[[j]]$Var2[i])
        int_crude_log_odds_matrix[,i] <- int[[j]]$beta[i] * xdata[,int_var1] * xdata[,int_var2]
      }
      
      int_crude_log_odds <- rowSums(int_crude_log_odds_matrix)
      
      # Linear combinations of effects
      
      crude_log_odds <- main_crude_log_odds + int_crude_log_odds
      
      # Identifying a relevant range of scaling factors (logistic distribution)
      s_max <- max(abs(crude_log_odds)) / log(0.51 / 0.49) # scale that gives all probabilities between 0.49 and 0.51 (expected c stat close to 0.5)
      s_min <- min(abs(crude_log_odds)) / log(0.99 / 0.01) # scale that gives all probabilities above 0.99 or below 0.01 (expected c stat close to 1)
      
      # Finding scaling factor that gives desired AUC (interval required to ensure optimisation works)
      argmax_scaling_factor <- stats::optimise(
        f = TuneCStatisticLogit,
        crude_log_odds = crude_log_odds,
        auc = ev_xy[j],
        lower = 1 / s_max, upper = 1 / s_min
      )
      scaling_factor <- argmax_scaling_factor$minimum
      
      # Applying scaling factor
      beta[, j] <- beta[, j] * scaling_factor
      log_odds <- crude_log_odds * scaling_factor
      
      # Calculating probabilities from log-odds (inverse logit)
      proba <- 1 / (1 + exp(-log_odds))
      
      # Sampling from Bernouilli distribution
      ydata[, j] <- stats::rbinom(n = n, size = 1, prob = proba)
    }
  }
  
  # Defining row and column names of ydata
  rownames(ydata) <- rownames(xdata)
  colnames(ydata) <- paste0("outcome", 1:q)
  
  # Defining column names of beta and theta
  colnames(beta) <- colnames(theta) <- colnames(ydata)
  
  # appending theta with theta_int
  
  theta_int <- matrix(0,nrow = ncol(xdata)*(ncol(xdata)-1)/2, ncol = q) # empty matrix
  rownames(theta_int) <- paste0(combn(colnames(xdata),2)[1,], sep=":", combn(colnames(xdata),2)[2,]) # interaction effect names
  
  for (i in 1:q){
    theta_int[,i][rownames(theta_int) %in% paste0(int[[i]]["Var1"][int[[i]]["theta"] == 1], sep=":", int[[i]][["Var2"]][int[[i]]["theta"] == 1])] <- int[[i]][["theta"]][int[[i]]["theta"] == 1] # find true interaction effects and replace with 1 
  }
  
  theta <- rbind(theta, theta_int)
  
  # appending beta with beta_int
  
  beta_int <- matrix(0,nrow = ncol(xdata)*(ncol(xdata)-1)/2, ncol = q) # empty matrix
  rownames(beta_int) <- paste0(combn(colnames(xdata),2)[1,], sep=":", combn(colnames(xdata),2)[2,]) # interaction effect names
  
  for (i in 1:q){
    beta_int[,i][rownames(beta_int) %in% paste0(int[[i]][["Var1"]][int[[i]]["beta"] != 0], sep=":", int[[i]][["Var2"]][int[[i]]["beta"] != 0])] <-  int[[i]]["beta"][int[[i]]["beta"] != 0] # find true interaction effects and replace with 1 
  }
  
  beta <- rbind(beta, beta_int)
  
  # Preparing parameters for output
  
  params <- c(pk, ev_xy)
  
  # Preparing the output
  out <- list(xdata = xdata, ydata = ydata, theta = theta, beta = beta)
  
  # Defining the class
  class(out) <- "simulation_interaction"
  
  return(out)
 }
 
### GLINTERNET implementation ----

PenalisedInteraction <- function(xdata,
                                 ydata,
                                 Lambda = NULL,
                                 family,
                                 ...){
  
  # Retrieving names of categorical and continuous variables
  xdata <- as.data.frame(xdata)
  cat_names <- colnames(xdata[,sapply(xdata, is.factor)])
  level_names <- do.call(rbind, lapply(xdata[,cat_names], levels))
  cont_names <- colnames(xdata[,sapply(xdata, is.numeric)])
  
  # Ordering categorical and continuous variables
  xdata_cat <- xdata[,cat_names]
  xdata_cont <- xdata[,cont_names]
  xdata <- cbind(xdata_cat, xdata_cont)
  
  # Calculating number of levels per variable
  numLevels <- apply(xdata,2, function(x) {
    if (is.factor(x))
      length(levels(x))
    else {1}})
  
  # Changing factor levels to fit the glinternet required format of categorical strings from 0:nlevels
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
      family = family
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
  index_main <- mapply(FUN = function(x,y){c(x, y)}, x = index_main_cat, y = index_main_cont)
  
  # Recording main effects based on indices
  selected_main <- as.data.frame(do.call(rbind,lapply(index_main, FUN = function(x){ifelse(attr(selected_main, "dimnames")[[2]] %in% x, 1,0)})))
  
  # Retrieving betas for main effects
  beta_main_cat <- lapply(mapply("[",mapply("[",coef)[2,])[1,], unlist)
  beta_main_cont <- lapply(mapply("[",mapply("[",coef)[2,])[2,], unlist)
  beta_main <- mapply(FUN = function(x,y){append(x, y)}, x = beta_main_cat, y = beta_main_cont)
  
  # Recording betas based on indices
  beta_main <- do.call(rbind,apply(as.matrix(Lambda),1, FUN = function(x){
    i <- which(Lambda == x)
    beta_main <- replace(selected_main[i,], which(selected_main[i,] == 1), beta_main[[i]])
    return(beta_main)
  })
  )
  
  # Generating matrix for Interaction effects
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
  
  # Retrieving indices matrix for interactions 
  index_int_catcat <- mapply("[",mapply("[",coef)[3,])[1,]
  index_int_contcont <- lapply(mapply("[",mapply("[",coef)[3,])[2,], FUN = function(x){x + length(xdata[,cat_names])})
  index_int_catcont <- lapply(mapply("[",mapply("[",coef)[3,])[3,], FUN = function(x){
    V1 <- x[,1]
    V2 <- x[,2] + length(xdata[,cat_names])
    out <- cbind(V1,V2)
    return(out)})
  index_int <- mapply(FUN = function(x,y,z){rbind(x,y,z)}, x = index_int_catcat, y = index_int_contcont, z = index_int_catcont)
  
  # Mapping glinternet indices to interaction indices
  index_int <- lapply(index_int, FUN = function(x){merge(x, index_int_dict, sort = F)})
  
  # Recording interaction effects based on indices
  selected_int<- as.data.frame(do.call(rbind,lapply(index_int, FUN = function(x){
    idx <- x[,3]
    selected_int <- ifelse(attr(matrix_interaction, "dimnames")[[2]] %in% idx, 1,0)
    return(selected_int)
    })))
  attr(selected_int, "assign") <- attr(matrix_interaction, "dimnames")[[2]] 
  
  # Retrieving betas for interaction effects
  beta_int_catcat <- lapply(mapply("[",mapply("[",coef)[4,])[1,], unlist)
  beta_int_contcont <- lapply(mapply("[",mapply("[",coef)[4,])[2,], unlist)
  beta_int_catcont <- lapply(mapply("[",mapply("[",coef)[4,])[3,], unlist)
  beta_int_full <- mapply(FUN = function(x,y,z){as.vector(c(x,y,z))}, x = beta_int_catcat, y = beta_int_contcont, z = beta_int_catcont)
  
  # Creating a dictionary for interaction betas
  beta_int_dict <- as.data.frame(attr(selected_int, "assign"))
  colnames(beta_int_dict) <- "V1:V2"
  beta_int_dict <- merge(index_int_dict, beta_int_dict)
  
  # Assigning indices to each list based on beta_dictinary and index_int
  beta_index_int <- lapply(index_int, FUN = function(x){merge(x, beta_int_dict, sort = F)})
  beta_int_map <- mapply(FUN = function(idx,beta){
    beta_int_map <- cbind(idx,beta)},
    idx = beta_index_int,
    beta = beta_int_full)
  beta_int_map <- lapply(beta_int_map, FUN = function(x){x[order(as.numeric(x[,3])),]})
  
  # Replacing based on column names
  beta_int <- matrix_interaction
  beta_int <- t(mapply(FUN = function(x,y){
    i <- which(Lambda == x)
    idx <- y[,3]
    beta <- y$beta
    beta_int <- replace(beta_int[i,], which(colnames(beta_int) %in% idx), beta)
  }, x = as.matrix(Lambda), y = beta_int_map))
  
  # Combining matrices
  selected <- cbind(selected_main, selected_int)
  beta_full <- cbind(beta_main, beta_int)
  
  # Changing row and column names
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
  
  # returning output 
  selected <- as.matrix(selected)
  beta_full <- as.matrix(beta_full)
  
  return(list(selected = selected, beta_full = beta_full))
  
}

### grid of penalty parameters using implementation glinternet ----
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

### cross-validated glinternet theta ----
glinternet.theta <- function(mymodel = NULL, xdata){
  
  # Retrieving names of categorical and continuous variables
  xdata <- as.data.frame(xdata)
  cat_names <- colnames(xdata[,sapply(xdata, is.factor)])
  level_names <- do.call(rbind, lapply(xdata[,cat_names], levels))
  cont_names <- colnames(xdata[,sapply(xdata, is.numeric)])
  
  # Ordering categorical and continuous variables
  xdata_cat <- xdata[,cat_names]
  xdata_cont <- xdata[,cont_names]
  xdata <- cbind(xdata_cat, xdata_cont)
  
  # Calculating number of levels per variable
  numLevels <- apply(xdata,2, function(x) {
    if (is.factor(x))
      length(levels(x))
    else {1}})
  
  # Retireving model coefficients
  
  coef <- coef(mymodel)
  
  # Retrieving Lambda 
  Lambda <- mymodel$lambda
  
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
  
  # Generating matrix for Interaction effects
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
  
  # Retrieving indices matrix for interactions 
  index_int_catcat <- mapply("[",mapply("[",coef)[3,])[1,]
  index_int_contcont <- lapply(mapply("[",mapply("[",coef)[3,])[2,], FUN = function(x){x + length(xdata[,cat_names])})
  index_int_catcont <- lapply(mapply("[",mapply("[",coef)[3,])[3,], FUN = function(x){
    V1 <- x[,1]
    V2 <- x[,2] + length(xdata[,cat_names])
    out <- cbind(V1,V2)
    return(out)})
  
  index_int <- mapply(FUN = function(x,y,z){rbind(x,y,z)}, x = index_int_catcat, y = index_int_contcont, z = index_int_catcont)
  
  # Mapping glinternet indices to interaction indices
  index_int <- lapply(index_int, FUN = function(x){merge(x, index_int_dict, sort = F)})
  
  # Recording interaction effects based on indices
  selected_int<- as.data.frame(do.call(rbind,lapply(index_int, FUN = function(x){
    idx <- x[,3]
    selected_int <- ifelse(attr(matrix_interaction, "dimnames")[[2]] %in% idx, 1,0)
    return(selected_int)
  })))
  
  attr(selected_int, "assign") <- attr(matrix_interaction, "dimnames")[[2]] 
  
  # Combining matrices
  selected <- cbind(selected_main, selected_int)
  # beta_full <- cbind(beta_main, beta_int)
  
  # Changing row and column names
  # colnames(beta_full) <- 
    colnames(selected) <- colnames(model.matrix( ~ 0 + .^2,
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
  # rownames(beta_full) <- paste0("s",0:(length(Lambda)-1))
  
  # returning output 
  selected <- as.matrix(selected)
  # beta_full <- as.matrix(beta_full)
  
  return(selected[2,])
  
}
 
glinternet.beta <- function(mymodel = NULL, xdata){
  
  # Retrieving names of categorical and continuous variables
  xdata <- as.data.frame(xdata)
  cat_names <- colnames(xdata[,sapply(xdata, is.factor)])
  level_names <- do.call(rbind, lapply(xdata[,cat_names], levels))
  cont_names <- colnames(xdata[,sapply(xdata, is.numeric)])
  
  # Ordering categorical and continuous variables
  xdata_cat <- xdata[,cat_names]
  xdata_cont <- xdata[,cont_names]
  xdata <- cbind(xdata_cat, xdata_cont)
  
  # Calculating number of levels per variable
  numLevels <- apply(xdata,2, function(x) {
    if (is.factor(x))
      length(levels(x))
    else {1}})
  
  # Retireving model coefficients
  coef <- coef(mymodel)
  
  # Retrieving Lambda 
  Lambda <- mymodel$lambda
  
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
  index_main <- mapply(FUN = function(x,y){c(x, y)}, x = index_main_cat, y = index_main_cont)
  
  # Recording main effects based on indices
  selected_main <- as.data.frame(do.call(rbind,lapply(index_main, FUN = function(x){ifelse(attr(selected_main, "dimnames")[[2]] %in% x, 1,0)})))
  
  # Retrieving betas for main effects
  beta_main_cat <- lapply(mapply("[",mapply("[",coef)[2,])[1,], unlist)
  beta_main_cont <- lapply(mapply("[",mapply("[",coef)[2,])[2,], unlist)
  beta_main <- mapply(FUN = function(x,y){append(x, y)}, x = beta_main_cat, y = beta_main_cont)
  
  # Recording betas based on indices
  beta_main <- do.call(rbind,apply(as.matrix(Lambda),1, FUN = function(x){
    i <- which(Lambda == x)
    beta_main <- replace(selected_main[i,], which(selected_main[i,] == 1), beta_main[[i]])
    return(beta_main)
  })
  )
  
  # Generating matrix for Interaction effects
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
  
  # Retrieving indices matrix for interactions 
  index_int_catcat <- mapply("[",mapply("[",coef)[3,])[1,]
  index_int_contcont <- lapply(mapply("[",mapply("[",coef)[3,])[2,], FUN = function(x){x + length(xdata[,cat_names])})
  index_int_catcont <- lapply(mapply("[",mapply("[",coef)[3,])[3,], FUN = function(x){
    V1 <- x[,1]
    V2 <- x[,2] + length(xdata[,cat_names])
    out <- cbind(V1,V2)
    return(out)})
  index_int <- mapply(FUN = function(x,y,z){rbind(x,y,z)}, x = index_int_catcat, y = index_int_contcont, z = index_int_catcont)
  
  # Mapping glinternet indices to interaction indices
  index_int <- lapply(index_int, FUN = function(x){merge(x, index_int_dict, sort = F)})
  
  # Recording interaction effects based on indices
  selected_int<- as.data.frame(do.call(rbind,lapply(index_int, FUN = function(x){
    idx <- x[,3]
    selected_int <- ifelse(attr(matrix_interaction, "dimnames")[[2]] %in% idx, 1,0)
    return(selected_int)
  })))
  attr(selected_int, "assign") <- attr(matrix_interaction, "dimnames")[[2]] 
  
  # Retrieving betas for interaction effects
  beta_int_catcat <- lapply(mapply("[",mapply("[",coef)[4,])[1,], unlist)
  beta_int_contcont <- lapply(mapply("[",mapply("[",coef)[4,])[2,], unlist)
  beta_int_catcont <- lapply(mapply("[",mapply("[",coef)[4,])[3,], unlist)
  beta_int_full <- mapply(FUN = function(x,y,z){as.vector(c(x,y,z))}, x = beta_int_catcat, y = beta_int_contcont, z = beta_int_catcont)
  
  # Creating a dictionary for interaction betas
  beta_int_dict <- as.data.frame(attr(selected_int, "assign"))
  colnames(beta_int_dict) <- "V1:V2"
  beta_int_dict <- merge(index_int_dict, beta_int_dict)
  
  # Assigning indices to each list based on beta_dictinary and index_int
  beta_index_int <- lapply(index_int, FUN = function(x){merge(x, beta_int_dict, sort = F)})
  beta_int_map <- mapply(FUN = function(idx,beta){
    beta_int_map <- cbind(idx,beta)},
    idx = beta_index_int,
    beta = beta_int_full)
  beta_int_map <- lapply(beta_int_map, FUN = function(x){x[order(as.numeric(x[,3])),]})
  
  # Replacing based on column names
  beta_int <- matrix_interaction
  beta_int <- t(mapply(FUN = function(x,y){
    i <- which(Lambda == x)
    idx <- y[,3]
    beta <- y$beta
    beta_int <- replace(beta_int[i,], which(colnames(beta_int) %in% idx), beta)
  }, x = as.matrix(Lambda), y = beta_int_map))
  
  # Combining matrices
  selected <- cbind(selected_main, selected_int)
  beta_full <- cbind(beta_main, beta_int)
  
  # Changing row and column names
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
  
  # returning output 
  selected <- as.matrix(selected)
  beta_full <- as.matrix(beta_full)
  
  return(beta = beta_full)
  
}

### [fake] functions ----
  
  MatchingArguments <- function(extra_args, FUN) {
    if ("..." %in% names(formals(FUN))) {
      out <- extra_args
    } else {
      ids <- which(names(extra_args) %in% names(formals(FUN)))
      out <- extra_args[ids]
    }
    return(out)
  }
  
  TuneCStatisticLogit <- function(scaling_factor, crude_log_odds, auc = NULL) {
    log_odds <- crude_log_odds * scaling_factor
    proba <- 1 / (1 + exp(-log_odds))
    if (is.null(auc)) {
      out <- ExpectedConcordance(probabilities = proba)
    } else {
      tmpauc <- ExpectedConcordance(probabilities = proba)
      out <- abs(auc - tmpauc)
    }
    return(out)
  }
  
  TuneExplainedVarianceCor <- function(u, ev_xx = NULL, omega) {
    diag(omega) <- diag(omega) + u
    mycor <- stats::cov2cor(solve(omega))
    tmp_ev <- norm(mycor, type = "2") / ncol(mycor)
    if (is.null(ev_xx)) {
      out <- tmp_ev
    } else {
      out <- abs(tmp_ev - ev_xx)
    }
    return(out)
  }
  
  TuneExplainedVarianceCov <- function(u, ev_xx = NULL, lambda) {
    lambda <- lambda + u
    lambda_inv <- 1 / lambda
    tmp_ev <- max(lambda_inv) / sum(lambda_inv)
    if (is.null(ev_xx)) {
      out <- tmp_ev
    } else {
      out <- abs(tmp_ev - ev_xx)
    }
    return(out)
  }
  
  MakePositiveDefinite <- function(omega, pd_strategy = "diagonally_dominant",
                                   ev_xx = NULL, scale = TRUE, u_list = c(1e-10, 1),
                                   tol = .Machine$double.eps^0.25) {
    # Making positive definite using diagonally dominance
    if (pd_strategy == "diagonally_dominant") {
      # Constructing the diagonal as the sum of entries
      diag(omega) <- apply(abs(omega), 1, sum)
      lambda <- eigen(omega, only.values = TRUE)$values
    }
    # Making positive definite using eigendecomposition
    if (pd_strategy == "min_eigenvalue") {
      # Extracting smallest eigenvalue of omega_tilde
      lambda <- eigen(omega, only.values = TRUE)$values
      lambda0 <- abs(min(lambda))
      
      # Making the precision matrix positive semidefinite
      lambda <- lambda + lambda0
      diag(omega) <- lambda0
    }
    
    if (is.null(ev_xx)) {
      # Finding u that maximises the contrast
      if (min(u_list) != max(u_list)) {
        argmax_u <- stats::optimise(MaxContrast,
                                    omega = omega, maximum = TRUE, tol = tol,
                                    lower = min(u_list), upper = max(u_list)
        )
        u <- argmax_u$maximum
      } else {
        u <- min(u_list)
      }
    } else {
      # Finding extreme values
      if (scale) {
        max_ev <- TuneExplainedVarianceCor(u = min(u_list), omega = omega)
        min_ev <- TuneExplainedVarianceCor(u = max(u_list), omega = omega)
      } else {
        max_ev <- TuneExplainedVarianceCov(u = min(u_list), lambda = lambda)
        min_ev <- TuneExplainedVarianceCov(u = max(u_list), lambda = lambda)
      }
      
      # Finding u corresponding to the required proportion of explained variance
      if ((ev_xx <= min_ev) | (ev_xx >= max_ev)) {
        if (ev_xx <= min_ev) {
          u <- max(u_list)
          if (ev_xx < min_ev) {
            message(paste0("The smallest proportion of explained variance by PC1 that can be obtained is ", round(min_ev, digits = 2), "."))
          }
        } else {
          u <- min(u_list)
          if (ev_xx > max_ev) {
            message(paste0("The largest proportion of explained variance by PC1 that can be obtained is ", round(max_ev, digits = 2), "."))
          }
        }
      } else {
        if (scale) {
          # Minimising the difference between requested and possible ev
          if (min(u_list) != max(u_list)) {
            argmin_u <- stats::optimise(TuneExplainedVarianceCor,
                                        omega = omega, ev_xx = ev_xx,
                                        lower = min(u_list), upper = max(u_list), tol = tol
            )
            u <- argmin_u$minimum
          } else {
            u <- min(u_list)
          }
        } else {
          # Minimising the difference between requested and possible ev
          if (min(u_list) != max(u_list)) {
            argmin_u <- stats::optimise(TuneExplainedVarianceCov,
                                        lambda = lambda, ev_xx = ev_xx,
                                        lower = min(u_list), upper = max(u_list), tol = tol
            )
            u <- argmin_u$minimum
          } else {
            u <- min(u_list)
          }
        }
      }
    }
    
    # Constructing the diagonal
    diag(omega) <- diag(omega) + u
    
    return(list(omega = omega, u = u))
  }
  
  MinWithinProba <- function(pk, nu_between = 0, nu_mat = NULL) {
    if (is.null(nu_mat)) {
      nu_mat <- diag(length(pk))
      nu_mat[upper.tri(nu_mat)] <- nu_between
      nu_mat[lower.tri(nu_mat)] <- nu_between
      diag(nu_mat) <- NA
    } else {
      if ((ncol(nu_mat) != length(pk)) & (nrow(nu_mat) != length(pk))) {
        stop("Arguments 'pk' and 'nu_mat' are not compatible. They correspond to different numbers of communities. The number of rows and columns in 'nu_mat' must be equal to the length of the vector 'pk'.")
      }
    }
    
    nu_within_min <- rep(NA, length(pk))
    for (k in 1:length(pk)) {
      nu_within_min[k] <- 1 / (pk[k] - 1) * (sum(nu_mat[k, -k] * pk[-k]) + 1 / pk[k])
    }
    
    return(nu_within_min)
  }
  
  ExpectedCommunities <- function(pk, nu_within = 0.1, nu_between = 0, nu_mat = NULL) {
    if (is.null(nu_mat)) {
      nu_mat <- diag(length(pk))
      nu_mat <- nu_mat * nu_within
      nu_mat[upper.tri(nu_mat)] <- nu_between
      nu_mat[lower.tri(nu_mat)] <- nu_between
    } else {
      if ((ncol(nu_mat) != length(pk)) & (nrow(nu_mat) != length(pk))) {
        stop("Arguments 'pk' and 'nu_mat' are not compatible. They correspond to different numbers of communities. The number of rows and columns in 'nu_mat' must be equal to the length of the vector 'pk'.")
      }
    }
    
    # Calculating expected sums of within and between degrees for community k
    d_within <- d_between <- rep(NA, length(pk))
    for (k in 1:length(pk)) {
      d_within[k] <- pk[k] * nu_mat[k, k] * (pk[k] - 1) # sum of within degrees in community k
      d_between[k] <- pk[k] * sum(nu_mat[k, -k] * pk[-k]) # sum of between degrees in community k
    }
    weak_community <- ifelse(d_within > d_between, yes = 1, no = 0)
    
    # Calculating expected number of edges within and between communities
    mystructure <- BlockStructure(pk)
    emat <- matrix(NA, nrow = nrow(mystructure), ncol = ncol(mystructure))
    for (k in unique(mystructure[upper.tri(mystructure, diag = TRUE)])) {
      if (k %in% diag(mystructure)) {
        i <- which(diag(mystructure) == k)
        emat[which(mystructure == k)] <- nu_mat[i, i] * pk[i] * (pk[i] - 1) / 2
      } else {
        tmpids <- which(mystructure == k, arr.ind = TRUE)
        i <- tmpids[1]
        j <- tmpids[2]
        emat[which(mystructure == k)] <- nu_mat[i, j] * pk[i] * pk[j]
      }
    }
    L <- sum(emat[upper.tri(emat, diag = TRUE)]) # total number of edges
    
    # Calculating community modularity
    M_C <- rep(NA, length(pk))
    for (k in 1:length(pk)) {
      M_C[k] <- emat[k, k] / L - ((d_within[k] + d_between[k]) / (2 * L))^2 # modularity of community k
    }
    
    # Calculating overall modularity
    M <- sum(M_C)
    
    return(list(
      total_within_degree_c = d_within,
      total_between_degree_c = d_between,
      weak_community = weak_community,
      total_number_edges_c = emat,
      modularity = M
    ))
  }
  
  
  Contrast <- function(mat, digits = 3) {
    return(length(unique(round(as.vector(abs(mat)), digits = digits))))
  }
  
  MaxContrast <- function(u, omega, digits = 3) {
    diag(omega) <- diag(omega) + u
    return(Contrast(stats::cov2cor(solve(omega)), digits = digits))
  }
  
  SimulateSymmetricMatrix <- function(pk = 10,
                                      v_within = c(0.5, 1), v_between = c(0, 0.1),
                                      v_sign = c(-1, 1), continuous = FALSE) {
    # Creating matrix with block indices
    bigblocks <- BlockMatrix(pk)
    bigblocks_vect <- bigblocks[upper.tri(bigblocks)]
    
    # Making as factor to allow for groups with 1 variable (for clustering)
    bigblocks_vect <- factor(bigblocks_vect, levels = seq(1, max(bigblocks)))
    block_ids <- unique(as.vector(bigblocks))
    
    # Building absolute v matrix
    v <- bigblocks
    v_vect <- v[upper.tri(v)]
    for (k in block_ids) {
      if (k %in% v_vect) {
        if (k %in% unique(diag(bigblocks))) {
          if (continuous) {
            v_vect[bigblocks_vect == k] <- stats::runif(sum(bigblocks_vect == k), min = min(v_within), max = max(v_within))
          } else {
            v_vect[bigblocks_vect == k] <- base::sample(v_within, size = sum(bigblocks_vect == k), replace = TRUE)
          }
        } else {
          if (continuous) {
            v_vect[bigblocks_vect == k] <- stats::runif(sum(bigblocks_vect == k), min = min(v_between), max = max(v_between))
          } else {
            v_vect[bigblocks_vect == k] <- base::sample(v_between, size = sum(bigblocks_vect == k), replace = TRUE)
          }
        }
      }
    }
    
    # Sampling the sign of precision entries
    v_vect <- v_vect * base::sample(sort(unique(v_sign)), size = length(v_vect), replace = TRUE)
    
    # Building v matrix
    diag(v) <- 0
    v[upper.tri(v)] <- v_vect
    v[lower.tri(v)] <- 0
    v <- v + t(v)
    
    return(v)
  }
  
  BlockDiagonal <- function(pk) {
    bigblocks <- BlockMatrix(pk)
    bigblocks[!bigblocks %in% diag(bigblocks)] <- 0
    bigblocks[bigblocks %in% diag(bigblocks)] <- 1
    return(bigblocks)
  }
  
  BlockStructure <- function(pk) {
    blocks <- BlockMatrix(pk = rep(1, length(pk)))
    
    return(blocks)
  }
  
  BlockMatrix <- function(pk) {
    nblocks <- sum(upper.tri(matrix(NA, ncol = length(pk), nrow = length(pk)), diag = TRUE))
    blocks <- matrix(NA, nrow = length(pk), ncol = length(pk))
    blocks[upper.tri(blocks, diag = TRUE)] <- 1:nblocks
    
    mybreaks <- c(0, cumsum(pk))
    bigblocks <- matrix(ncol = sum(pk), nrow = sum(pk))
    row_id_start <- matrix(mybreaks[row(blocks)], ncol = length(pk)) + 1
    row_id_end <- matrix(mybreaks[row(blocks) + 1], ncol = length(pk))
    col_id_start <- matrix(mybreaks[col(blocks)], ncol = length(pk)) + 1
    col_id_end <- matrix(mybreaks[col(blocks) + 1], ncol = length(pk))
    
    row_id_start <- row_id_start[upper.tri(row_id_start, diag = TRUE)]
    row_id_end <- row_id_end[upper.tri(row_id_end, diag = TRUE)]
    col_id_start <- col_id_start[upper.tri(col_id_start, diag = TRUE)]
    col_id_end <- col_id_end[upper.tri(col_id_end, diag = TRUE)]
    
    for (block_id in blocks[upper.tri(blocks, diag = TRUE)]) {
      ids <- rbind(
        expand.grid(
          row_id_start[block_id]:row_id_end[block_id],
          col_id_start[block_id]:col_id_end[block_id]
        ),
        expand.grid(
          col_id_start[block_id]:col_id_end[block_id],
          row_id_start[block_id]:row_id_end[block_id]
        )
      )
      bigblocks[as.matrix(ids)] <- block_id
    }
    
    return(bigblocks)
  }
  
  SimulatePrecision <- function(pk = NULL, theta,
                                v_within = c(0.5, 1), v_between = c(0, 0.1),
                                v_sign = c(-1, 1), continuous = TRUE,
                                pd_strategy = "diagonally_dominant", ev_xx = NULL, scale = TRUE,
                                u_list = c(1e-10, 1), tol = .Machine$double.eps^0.25) {
    # Checking inputs and defining pk
    if (is.null(pk)) {
      pk <- ncol(theta)
    } else {
      if (sum(pk) != ncol(theta)) {
        stop("Arguments 'pk' and 'theta' are not consistent. The sum of 'pk' entries must be equal to the number of rows and columns in 'theta'.")
      }
    }
    
    # Checking the choice of pd_strategy
    if (!pd_strategy %in% c("diagonally_dominant", "min_eigenvalue")) {
      stop("Invalid input for argument 'pd_strategy'. Possible values are: 'diagonally_dominant' or 'min_eigenvalue'.")
    }
    
    # Checking other input values
    if (any((v_within < 0) | (v_within > 1))) {
      stop("Invalid input for argument 'v_within'. Values must be between 0 and 1.")
    }
    if (any((v_between < 0) | (v_between > 1))) {
      stop("Invalid input for argument 'v_between'. Values must be between 0 and 1.")
    }
    if (any(!v_sign %in% c(-1, 1))) {
      stop("Invalid input for argument 'v_sign'. Possible values are -1 and 1.")
    }
    
    # Ensuring that v values are lower than or equal to 1
    if (any(abs(v_within) > 1)) {
      v_within <- v_within / max(abs(v_within))
      message("The values provided in 'v_within' have been re-scaled to be lower than or equal to 1 in absolute value.")
    }
    
    # Ensuring that diagonal entries of theta are zero
    diag(theta) <- 0
    
    # Building v matrix
    v <- SimulateSymmetricMatrix(
      pk = pk, v_within = v_within, v_between = v_between,
      v_sign = v_sign, continuous = continuous
    )
    
    # Filling off-diagonal entries of the precision matrix
    omega_tilde <- theta * v
    
    # Ensuring positive definiteness
    omega_pd <- MakePositiveDefinite(
      omega = omega_tilde, pd_strategy = pd_strategy,
      ev_xx = ev_xx, scale = scale, u_list = u_list, tol = tol
    )
    
    # Returning the output
    return(omega_pd)
  }
  
SimulateAdjacency <- function(pk = 10,
                              implementation = HugeAdjacency,
                              topology = "random",
                              nu_within = 0.1,
                              nu_between = 0,
                              nu_mat = NULL,
                              ...) {
    # Storing all arguments
    args <- c(mget(ls()), list(...))
    
    # Checking the inputs
    if (topology != "random") {
      if (length(pk) > 1) {
        pk <- sum(pk)
        warning(paste0("Multi-block simulations are only allowed with topology='random'. Argument 'pk' has been set to ", pk, "."))
      }
    }
    
    # Creating the matrix of probabilities
    if (is.null(nu_mat)) {
      nu_mat <- diag(length(pk)) * nu_within
      nu_mat[upper.tri(nu_mat)] <- nu_between
      nu_mat[lower.tri(nu_mat)] <- nu_between
    } else {
      if ((ncol(nu_mat) != length(pk)) & (nrow(nu_mat) != length(pk))) {
        stop("Arguments 'pk' and 'nu_mat' are not compatible. They correspond to different numbers of communities. The number of rows and columns in 'nu_mat' must be equal to the length of the vector 'pk'.")
      }
    }
    
    # Creating matrix with block indices
    bigblocks <- BlockMatrix(pk)
    bigblocks_vect <- bigblocks[upper.tri(bigblocks)]
    
    # Making as factor to allow for groups with 1 variable (for clustering)
    bigblocks_vect <- factor(bigblocks_vect, levels = seq(1, max(bigblocks)))
    block_ids <- unique(as.vector(bigblocks))
    
    # Creating matrix with block structure
    blocks <- BlockStructure(pk)
    
    # Identifying relevant arguments
    if (!"..." %in% names(formals(implementation))) {
      ids <- which(names(args) %in% names(formals(implementation)))
      args <- args[ids]
    }
    
    # Simulation of the adjacency matrix
    if ("nu" %in% names(formals(implementation))) {
      if (length(pk) > 1) {
        # Initialising theta
        theta <- matrix(0, nrow = sum(pk), ncol = sum(pk))
        theta_vect <- theta[upper.tri(theta)]
        
        # # Allowing for different densities in within and between blocks
        # theta_w <- do.call(implementation, args = c(args, list(nu = nu_within)))
        # theta_w_vect <- theta_w[upper.tri(theta_w)]
        # theta_b <- do.call(implementation, args = c(args, list(nu = nu_between)))
        # theta_b_vect <- theta_b[upper.tri(theta_b)]
        
        # Filling within and between blocks
        for (k in block_ids) {
          tmpids <- which(blocks == k, arr.ind = TRUE)
          i <- tmpids[1]
          j <- tmpids[2]
          theta_w <- do.call(implementation, args = c(args, list(nu = nu_mat[i, j])))
          theta_w_vect <- theta_w[upper.tri(theta_w)]
          theta_vect[bigblocks_vect == k] <- theta_w_vect[bigblocks_vect == k]
          # if (k %in% unique(diag(bigblocks))) {
          #   theta_vect[bigblocks_vect == k] <- theta_w_vect[bigblocks_vect == k]
          # } else {
          #   theta_vect[bigblocks_vect == k] <- theta_b_vect[bigblocks_vect == k]
          # }
        }
        theta[upper.tri(theta)] <- theta_vect
        theta <- theta + t(theta)
      } else {
        theta <- do.call(implementation, args = c(args, list(nu = nu_within)))
      }
    } else {
      theta <- do.call(implementation, args = c(args))
    }
    
    # Ensuring the adjacency matrix is symmetric (undirected graph) with no self-loops
    theta <- ifelse(theta + t(theta) != 0, yes = 1, no = 0)
    diag(theta) <- 0
    
    # Setting variable names
    colnames(theta) <- rownames(theta) <- paste0("var", 1:ncol(theta))
    
    # Defining the class
    class(theta) <- c("matrix", "adjacency_matrix")
    
    return(theta)
  }
  
SimulateGraphical <- function(n = 100, pk = 10, theta = NULL,
                              implementation = HugeAdjacency, topology = "random",
                              nu_within = 0.1, nu_between = NULL, nu_mat = NULL,
                              v_within = c(0.5, 1), v_between = c(0.1, 0.2),
                              v_sign = c(-1, 1), continuous = TRUE,
                              pd_strategy = "diagonally_dominant", ev_xx = NULL, scale_ev = TRUE,
                              u_list = c(1e-10, 1), tol = .Machine$double.eps^0.25,
                              scale = TRUE, output_matrices = FALSE, ...) {
    # Defining number of nodes
    p <- sum(pk)
    if (!is.null(theta)) {
      if (ncol(theta) != p) {
        p <- pk <- ncol(theta)
      }
    }
    
    # Defining the between-block density
    if (is.null(nu_between)) {
      nu_between <- nu_within
    }
    
    # Building adjacency matrix
    if (is.null(theta)) {
      theta <- SimulateAdjacency(
        pk = pk,
        implementation = implementation, topology = topology,
        nu_within = nu_within, nu_between = nu_between, nu_mat = nu_mat, ...
      )
    }
    
    # Simulation of a precision matrix
    out <- SimulatePrecision(
      pk = pk, theta = theta,
      v_within = v_within, v_between = v_between,
      v_sign = v_sign, continuous = continuous,
      pd_strategy = pd_strategy, ev_xx = ev_xx, scale = scale_ev,
      u_list = u_list, tol = tol
    )
    omega <- out$omega
    
    # Computing the covariance matrix
    if (scale) {
      sigma <- stats::cov2cor(solve(omega))
    } else {
      sigma <- solve(omega)
    }
    
    # Computing the partial correlation matrix
    if (output_matrices) {
      phi <- -stats::cov2cor(omega) + 2 * diag(ncol(omega))
    }
    
    # Simulating data from multivariate normal distribution
    x <- MASS::mvrnorm(n, rep(0, p), sigma)
    colnames(x) <- paste0("var", 1:ncol(x))
    rownames(x) <- paste0("obs", 1:nrow(x))
    
    # Defining the class of theta
    class(theta) <- c("matrix", "adjacency_matrix")
    
    if (output_matrices) {
      out <- list(
        data = x, theta = theta,
        omega = omega, phi = phi, sigma = sigma,
        u = out$u
      )
    } else {
      out <- list(data = x, theta = theta)
    }
    
    # Defining the class
    class(out) <- "simulation_graphical_model"
    
    return(out)
  }  


HugeAdjacency <- function(pk = 10, topology = "random", nu = 0.1, ...) {
  # Storing extra arguments
  extra_args <- list(...)
  
  # Extracting relevant extra arguments
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = huge::huge.generator)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("n", "d", "prob", "graph", "verbose")]
  
  # Running simulation model
  mymodel <- do.call(huge::huge.generator, args = c(
    list(
      n = 2, d = sum(pk), prob = nu,
      graph = topology, verbose = FALSE
    ),
    tmp_extra_args
  ))
  theta <- as.matrix(mymodel$theta)
  
  # Re-organising the variables to avoid having centrality related to variable ID (e.g. for scale-free models)
  ids <- sample(ncol(theta))
  theta <- theta[ids, ids]
  
  return(theta)
}


SamplePredictors <- function(pk, q = NULL, nu = 0.1, orthogonal = TRUE) {
  # Definition of the number of outcome variables
  if (is.null(q)) {
    q <- length(pk)
  }
  if (length(nu) != q) {
    nu <- rep(nu[1], q)
  }
  
  # Simulation of the binary status for true predictors
  theta <- matrix(0, nrow = q, ncol = sum(pk))
  for (k in 1:q) {
    if (orthogonal) {
      if (k > 1) {
        ids <- seq(cumsum(pk)[k - 1] + 1, cumsum(pk)[k])
      } else {
        ids <- seq(1, cumsum(pk)[k])
      }
      theta[k, ids] <- stats::rbinom(pk[k], size = 1, prob = nu[k])
      
      # Introducing at least one true predictor
      if (sum(theta[k, ids]) == 0) {
        theta[k, sample(ids, size = 1)] <- 1
      }
    } else {
      theta[k, ] <- stats::rbinom(sum(pk), size = 1, prob = nu[k])
      theta[k, k] <- 1
    }
  }
  
  return(t(theta))
}


### [sharp] functions ----
Refit <- function(xdata, ydata, stability = NULL,
                  family = NULL, implementation = NULL,
                  Lambda = NULL, seed = 1,
                  verbose = TRUE, ...) {
  # Storing extra arguments
  extra_args <- list(...)
  if ("check_input" %in% names(extra_args)) {
    check_input <- extra_args$check_input
    extra_args <- extra_args[!names(extra_args) %in% "check_input"]
  } else {
    check_input <- TRUE
  }
  
  # Defining the type of model (PLS vs regression)
  use_pls <- FALSE
  
  # Checking input
  if (!is.null(stability)) {
    if (!inherits(stability, c("variable_selection", "bi_selection"))) {
      stop("Argument 'stability' is not of class 'variable_selection' or 'bi_selection'. This function can only be applied on the output of (i) VariableSelection() or (ii) BiSelection() for PLS models.")
    }
    if (inherits(stability, "bi_selection")) {
      # Checking mixOmics package is installed
      CheckPackageInstalled("mixOmics")
      use_pls <- TRUE
      
      if (!stability$methods$family %in% c("gaussian")) {
        stop("This function can only be applied with the 'gaussian' family for PLS models.")
      }
    } else {
      if (!stability$methods$family %in% c("gaussian", "cox", "binomial", "multinomial")) {
        stop("This function can only be applied with the following families for regression models: 'gaussian', 'cox', 'binomial', or 'multinomial'.")
      }
    }
    if (!is.null(family)) {
      if (family != stability$methods$family) {
        warning(paste0("Arguments 'stability' and 'family' are not consistent. The family specified in argument stability was used: ", stability$methods$family))
      }
    }
    family <- stability$methods$family
  } else {
    if (is.null(family)) {
      stop("Argument 'family' must be provided. Possible values are: 'gaussian', 'cox', 'binomial', or 'multinomial'.")
    }
  }
  
  # Object preparation, error and warning messages
  if (check_input) {
    CheckDataRegression(
      xdata = xdata, ydata = ydata, family = family, verbose = verbose
    )
  }
  
  if (use_pls) {
    # Refitting the PLS model
    mymodel <- PLS(
      xdata = xdata, ydata = ydata,
      selectedX = stability$selectedX,
      selectedY = stability$selectedY,
      family = family, ncomp = NULL,
      scale = stability$methods$scale
    )
  } else {
    # Extracting the stably selected predictors
    if (is.null(stability)) {
      selected <- rep(1, ncol(xdata))
      names(selected) <- colnames(xdata)
    } else {
      selected <- SelectedVariables(stability)
    }
    
    # Defining predictors for the model (including un-penalised)
    ids <- c(
      names(selected)[which(selected == 1)],
      colnames(xdata)[!colnames(xdata) %in% names(selected)]
    )
    xdata <- xdata[, ids, drop = FALSE]
    
    if (is.null(implementation)) {
      # Writing model formula
      ids <- gsub("`", "", ids)
      colnames(xdata) <- gsub("`", "", colnames(xdata))
      if (length(ids) == 0) {
        message("No stably selected variables. Running a model with intercept only.")
        myformula <- stats::as.formula("ydata ~ 1")
      } else {
        myformula <- stats::as.formula(paste0("ydata ~ ", paste(paste0("`", ids, "`"), collapse = " + ")))
      }
      
      # Defining penalisation
      penalised <- TRUE
      if (ncol(xdata) == 1) {
        penalised <- FALSE
      }
      
      # Preparing the model
      if (penalised) {
        # Extracting relevant extra arguments
        tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = glmnet::cv.glmnet)
        tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("x", "y", "lambda", "alpha", "family", "type.multinomial")]
        
        # Setting seed for reproducibility
        withr::local_seed(seed)
        
        # Running model
        if ((is.null(Lambda)) | (length(Lambda) > 1)) {
          mymodel <- do.call(glmnet::cv.glmnet, args = c(
            list(
              x = xdata,
              y = ydata,
              family = family,
              alpha = 0,
              lambda = Lambda,
              type.multinomial = "grouped"
            ),
            tmp_extra_args
          ))
        } else {
          mymodel <- do.call(glmnet::glmnet, args = c(
            list(
              x = xdata,
              y = ydata,
              family = family,
              alpha = 0,
              lambda = Lambda,
              type.multinomial = "grouped"
            ),
            tmp_extra_args
          ))
        }
      } else {
        # Recalibration for linear regression
        if (family == "gaussian") {
          tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = stats::lm)
          tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("formula", "data")]
          mymodel <- do.call(stats::lm, args = c(
            list(
              formula = myformula,
              data = as.data.frame(xdata)
            ),
            tmp_extra_args
          ))
        }
        
        # Recalibration for Cox regression
        if (family == "cox") {
          ids1 <- which(names(extra_args) %in% names(formals(survival::coxph)))
          ids2 <- which(names(extra_args) %in% names(formals(survival::coxph.control)))
          tmp_extra_args <- extra_args[unique(c(ids1, ids2))]
          tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("formula", "data")]
          ydata <- survival::Surv(time = ydata[, 1], event = ydata[, 2])
          mymodel <- do.call(survival::coxph, args = c(
            list(
              formula = myformula,
              data = as.data.frame(xdata)
            ),
            tmp_extra_args
          ))
        }
        
        # Recalibration for logistic regression
        if (family == "binomial") {
          ids1 <- which(names(extra_args) %in% names(formals(stats::glm)))
          ids2 <- which(names(extra_args) %in% names(formals(stats::glm.control)))
          tmp_extra_args <- extra_args[unique(c(ids1, ids2))]
          tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("formula", "data", "family")]
          suppressWarnings({
            mymodel <- do.call(stats::glm, args = c(
              list(
                formula = myformula,
                data = as.data.frame(xdata),
                family = stats::binomial(link = "logit")
              ),
              tmp_extra_args
            ))
          })
        }
        
        # Recalibration for multinomial regression
        if (family == "multinomial") {
          tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = nnet::multinom)
          tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("formula", "data", "trace")]
          mymodel <- do.call(nnet::multinom, args = c(
            list(
              formula = myformula,
              data = as.data.frame(xdata),
              trace = FALSE
            ),
            tmp_extra_args
          ))
        }
      }
    } else {
      tmp_extra_args <- extra_args[!names(extra_args) %in% c("xdata", "ydata", "family")]
      xdata <- xdata[, ids, drop = FALSE]
      mymodel <- do.call(implementation, args = c(
        list(
          xdata = xdata,
          ydata = ydata,
          family = family
        ),
        tmp_extra_args
      ))
    }
  }
  
  return(mymodel)
}

ExplanatoryPerformance <- function(xdata, ydata, new_xdata = NULL, new_ydata = NULL,
                                   stability = NULL, family = NULL,
                                   implementation = NULL, prediction = NULL, resampling = "subsampling",
                                   K = 1, tau = 0.8, seed = 1,
                                   n_thr = NULL,
                                   time = 1000,
                                   verbose = FALSE, ...) {
  # Checking the inputs
  if (!is.null(stability)) {
    if (!inherits(stability, "variable_selection")) {
      stop("Argument 'stability' is not of class 'variable_selection'. This function can only be applied on the output of VariableSelection().")
    }
    if (!stability$methods$family %in% c("cox", "binomial", "gaussian")) {
      stop("This function can only be applied with the following families: 'binomial', 'cox' or 'gaussian'.")
    }
    if (!is.null(family)) {
      if (family != stability$methods$family) {
        warning(paste0("Arguments 'stability' and 'family' are not consistent. The family specified in argument stability was used: ", stability$methods$family))
      }
    }
    family <- stability$methods$family
  } else {
    if (is.null(family)) {
      stop("Argument 'family' must be provided. Possible values are: 'gaussian', 'cox' or 'binomial'.")
    }
  }
  if (!is.null(new_xdata)) {
    if (is.null(new_ydata)) {
      stop("Argument 'new_ydata' must be provided if 'new_xdata' is provided.")
    }
    K <- 1
  }
  
  # Object preparation, error and warning messages
  CheckDataRegression(
    xdata = xdata, ydata = ydata, family = family, verbose = verbose
  )
  if (!is.null(new_xdata)) {
    xtrain <- xdata
    ytrain <- ydata
    CheckDataRegression(
      xdata = new_xdata, ydata = new_ydata, family = family, verbose = verbose
    )
    xtest <- xdata # variable updated by CheckDataRegression()
    ytest <- ydata
  }
  
  # Defining the metric to use
  if (family == "binomial") {
    metric <- "roc"
  }
  if (family == "cox") {
    metric <- "concordance"
  }
  if (family == "gaussian") {
    metric <- "q2"
  }
  
  # Setting seed for reproducibility
  withr::local_seed(seed)
  
  # Running the subsampling iterations
  iter <- 0
  for (k in 1:K) {
    iter <- iter + 1
    
    if (is.null(new_xdata)) {
      # Balanced training/test split
      ids_test <- Resample(data = ydata, tau = 1 - tau, family = family, resampling = resampling, ...)
      xtrain <- xdata[-ids_test, , drop = FALSE]
      ytrain <- ydata[-ids_test, , drop = FALSE]
      xtest <- xdata[ids_test, , drop = FALSE]
      ytest <- ydata[ids_test, , drop = FALSE]
    }
    
    # Recalibration from stability selection model
    refitted <- Refit(
      xdata = xtrain, ydata = ytrain,
      stability = stability,
      implementation = implementation,
      family = family,
      check_input = FALSE,
      ...
    )
    
    if (is.null(implementation)) {
      # Initialising matrix of beta coefficients
      if (iter == 1) {
        if (family %in% c("gaussian", "binomial", "cox")) {
          tmpbeta <- as.vector(stats::coef(refitted))
          Beta <- matrix(NA, nrow = K, ncol = length(tmpbeta))
          colnames(Beta) <- names(tmpbeta)
          rownames(Beta) <- paste0("iter", 1:K)
        }
      }
      
      # Storing beta coefficients
      if (family %in% c("gaussian", "binomial", "cox")) {
        Beta[iter, ] <- as.vector(stats::coef(refitted))
      }
      
      if (ncol(xtrain) == 1) {
        # Predictions from logistic models
        if (tolower(metric) == "roc") {
          suppressWarnings({
            predicted <- stats::predict.glm(refitted, newdata = as.data.frame(xtest))
          })
        }
        
        # Predictions from linear models
        if (tolower(metric) == "q2") {
          predicted <- stats::predict.lm(refitted, newdata = as.data.frame(xtest))
        }
      } else {
        ids_predictors <- rownames(stats::coef(refitted))
        ids_predictors <- ids_predictors[which(ids_predictors %in% colnames(xtest))]
        predicted <- stats::predict(refitted, newx = as.matrix(xtest[, ids_predictors]))
      }
    } else {
      if (is.null(prediction)) {
        stop("Argument 'prediction' has to be provided if 'implementation' is provided. It must be a function that takes the output of 'implementation' as argument.")
      }
      predicted <- do.call(prediction, args = list(xdata = xtest, model = refitted))
    }
    
    # Performing ROC analyses
    if (tolower(metric) == "roc") {
      # ROC analysis
      roc <- ROC(predicted = predicted, observed = ytest, n_thr = length(predicted) - 1)
      
      # Initialisation of the object
      if (iter == 1) {
        FPR <- TPR <- matrix(NA, nrow = K, ncol = length(roc$TPR))
        AUC <- rep(NA, K)
      }
      
      # Storing the metrics
      FPR[iter, ] <- roc$FPR
      TPR[iter, ] <- roc$TPR
      AUC[iter] <- roc$AUC
    }
    
    # Performing Q-squared analyses
    if (tolower(metric) == "q2") {
      # Initialisation of the object
      if (iter == 1) {
        Q_squared <- rep(NA, K)
      }
      
      # Computing the Q-squared
      Q_squared[iter] <- stats::cor(predicted, ytest)^2
    }
    
    # Performing concordance analyses
    if (tolower(metric) == "concordance") {
      # Computing the concordance index for given times
      if (length(stats::coef(refitted)) == 1) {
        predicted <- stats::predict(refitted, newdata = as.data.frame(xtest), type = "lp")
      } else {
        predicted <- stats::predict(refitted, newx = as.matrix(xtest), type = "link")
      }
      ytest[which(ytest[, "time"] > time), "status"] <- 0
      cstat <- glmnet::Cindex(pred = predicted, y = ytest)
      if (iter == 1) {
        cindex <- rep(NA, K)
      }
      cindex[iter] <- cstat
    }
  }
  
  # Preparing the output
  if (tolower(metric) == "roc") {
    out <- list(FPR = FPR, TPR = TPR, AUC = AUC)
  }
  
  if (tolower(metric) == "concordance") {
    out <- list(concordance = cindex)
  }
  
  if (tolower(metric) == "q2") {
    out <- list(Q_squared = Q_squared)
  }
  
  if (is.null(implementation)) {
    out <- c(out, Beta = list(Beta))
  }
  
  # Defining class
  class(out) <- "roc_band"
  
  return(out)
}

Incremental <- function(xdata, ydata, new_xdata = NULL, new_ydata = NULL,
                        stability = NULL, family = NULL,
                        implementation = NULL, prediction = NULL, resampling = "subsampling",
                        n_predictors = NULL,
                        K = 100, tau = 0.8, seed = 1,
                        n_thr = NULL,
                        time = 1000,
                        verbose = TRUE, ...) {
  # Checking the inputs
  if (!is.null(stability)) {
    if (!inherits(stability, "variable_selection")) {
      stop("Argument 'stability' is not of class 'variable_selection'. This function can only be applied on the output of VariableSelection().")
    }
    if (!is.null(family)) {
      if (family != stability$methods$family) {
        warning(paste0("Arguments 'stability' and 'family' are not consistent. The family specified in argument stability was used: ", stability$methods$family))
      }
    }
    family <- stability$methods$family
    if (!family %in% c("cox", "binomial", "gaussian")) {
      stop("This function can only be applied with the following families: 'binomial', 'cox' or 'gaussian'.")
    }
  } else {
    if (is.null(family)) {
      stop("Argument 'family' must be provided. Possible values are: 'binomial', 'cox' or 'gaussian'.")
    }
  }
  
  # Defining the number of predictors
  if (is.null(n_predictors)) {
    if (!is.null(stability)) {
      # Stopping at the calibrated model
      n_predictors <- sum(SelectedVariables(stability))
      
      # Adding the variables that are forced in the model with penalty.factor
      n_predictors <- n_predictors + sum(!colnames(xdata) %in% names(SelectedVariables(stability)))
    } else {
      n_predictors <- ncol(xdata)
    }
    if (n_predictors == 0) {
      n_predictors <- 10
    }
  }
  
  # Defining the order of inclusion in the model
  if (is.null(stability)) {
    myorder <- colnames(xdata) # order of the columns is used
    n_predictors <- min(n_predictors, ncol(xdata))
  } else {
    # Identifying unpenalised variables
    unpenalised_vars <- colnames(xdata)[!colnames(xdata) %in% names(SelectedVariables(stability))]
    n_unpenalised <- length(unpenalised_vars)
    n_predictors <- min(n_predictors, ncol(xdata) - n_unpenalised)
    
    # Including variables by order of decreasing selection proportions
    myorder <- names(SelectionProportions(stability))[sort.list(SelectionProportions(stability), decreasing = TRUE)]
    myorder <- myorder[1:n_predictors]
    
    # Including the variables that are forced in the model first by order of columns in the data
    myorder <- c(unpenalised_vars, myorder)
  }
  
  # Initialisation of the objects
  Beta <- list()
  if (family == "binomial") {
    TPR <- FPR <- AUC <- list()
  }
  if (family == "cox") {
    concordance <- list()
  }
  if (family == "gaussian") {
    Q_squared <- list()
  }
  
  if (verbose) {
    pb <- utils::txtProgressBar(style = 3)
  }
  
  for (k in 1:length(myorder)) {
    perf <- ExplanatoryPerformance(
      xdata = xdata[, myorder[1:k], drop = FALSE],
      ydata = ydata,
      new_xdata = new_xdata,
      new_ydata = new_ydata,
      stability = NULL,
      family = family,
      implementation = implementation,
      prediction = prediction,
      resampling = resampling,
      K = K, tau = tau, seed = seed,
      n_thr = n_thr,
      time = time,
      verbose = FALSE,
      ...
    )
    if (family == "binomial") {
      FPR <- c(FPR, list(perf$FPR))
      TPR <- c(TPR, list(perf$TPR))
      AUC <- c(AUC, list(perf$AUC))
    }
    if (family == "cox") {
      concordance <- c(concordance, list(perf$concordance))
    }
    if (family == "gaussian") {
      Q_squared <- c(Q_squared, list(perf$Q_squared))
    }
    Beta <- c(Beta, list(perf$Beta))
    
    if (verbose) {
      utils::setTxtProgressBar(pb, k / length(myorder))
    }
  }
  
  # Preparing the output
  if (family == "binomial") {
    out <- list(FPR = FPR, TPR = TPR, AUC = AUC)
  }
  if (family == "cox") {
    out <- list(concordance = concordance)
  }
  if (family == "gaussian") {
    out <- list(Q_squared = Q_squared)
  }
  
  # Adding beta coefficients
  if (is.null(implementation)) {
    out <- c(out, list(Beta = Beta))
  }
  
  # Adding variable names
  mynames <- myorder
  out <- c(out, names = list(mynames))
  
  # Adding stable selection status
  if (!is.null(stability)) {
    mystable <- rep(0, length(myorder))
    names(mystable) <- myorder
    mystable[1:max(which(myorder %in% names(which(SelectedVariables(stability) == 1))))] <- 1
    out <- c(out, stable = list(mystable))
  }
  
  # Defining class
  class(out) <- "incremental"
  
  return(out)
}

Combine <- function(stability1, stability2, include_beta = TRUE) {
  if (!inherits(stability1, c("variable_selection", "structural_model", "graphical_model", "clustering"))) {
    stop("Invalid inputs. This function only applies to outputs from VariableSelection(), StructuralModel(), GraphicalModel() or Clustering().")
  }
  if (class(stability1) != class(stability2)) {
    stop("Arguments 'stability1' and 'stability2' are not compatible. They must be generated from the same function.")
  }
  if (!all(is.na(stability1$Lambda))) {
    if (any(stability1$Lambda != stability2$Lambda)) {
      stop("Arguments 'stability1' and 'stability2' are not compatible. They were constructed using different Lambdas.")
    }
  }
  if (any(do.call(c, stability1$methods) != do.call(c, stability2$methods))) {
    stop("Arguments 'stability1' and 'stability2' are not compatible. They were constructed using different methods.")
  }
  if (inherits(stability1, "clustering")) {
    if (any(do.call(c, stability1$params[c("pk", "n", "tau")]) != do.call(c, stability2$params[c("pk", "n", "tau")]))) {
      stop("Arguments 'stability1' and 'stability2' are not compatible. They were constructed using different methods.")
    }
  } else {
    if (any(do.call(c, stability1$params[c("pk", "n", "tau", "PFER_thr", "FDP_thr")]) != do.call(c, stability2$params[c("pk", "n", "tau", "PFER_thr", "FDP_thr")]))) {
      stop("Arguments 'stability1' and 'stability2' are not compatible. They were constructed using different methods.")
    }
  }
  if (stability1$params$seed == stability2$params$seed) {
    stop("Arguments 'stability1' and 'stability2' were obtained using the same seed.")
  }
  if (any(stability1$sign != stability2$sign)) {
    stop("Arguments 'stability1' and 'stability2' are not compatible. They were obtained from different datasets.")
  }
  
  # Identifying the type of model (variable/pairs)
  if (inherits(stability1, "graphical_model")) {
    graph <- TRUE
  } else {
    graph <- FALSE
  }
  
  # Extracting the parameters (some are NULL depending on the method)
  nc <- stability1$nc
  Lambda <- stability1$Lambda
  mysign <- stability1$sign
  pk <- stability1$params$pk
  pi_list <- stability1$params$pi_list
  n_cat <- stability1$params$n_cat
  Sequential_template <- stability1$params$Sequential_template
  PFER_method <- stability1$methods$PFER_method
  PFER_thr <- stability1$params$PFER_thr
  FDP_thr <- stability1$params$FDP_thr
  mymethods <- stability1$methods
  myparams <- stability1$params
  
  # Creating matrix with block indices (specific to graphical models)
  nblocks <- 1
  N_blocks <- stability1$params$pk
  names(N_blocks) <- 1
  if (stability1$methods$type == "graphical_model") { # to avoid memory issues in high dimensional variable selection
    bigblocks <- BlockMatrix(pk)
    bigblocks_vect <- bigblocks[upper.tri(bigblocks)]
    N_blocks <- unname(table(bigblocks_vect))
    blocks <- unique(as.vector(bigblocks_vect))
    names(N_blocks) <- blocks
    nblocks <- max(blocks)
  }
  
  # Preparing the PFER and FDP thresholds
  if (length(PFER_thr) == 1) {
    PFER_thr_blocks <- ceiling(prop.table(N_blocks) * PFER_thr)
  } else {
    if (length(PFER_thr) == nblocks) {
      PFER_thr_blocks <- PFER_thr
    }
  }
  if (length(FDP_thr) == 1) {
    FDP_thr_blocks <- rep(FDP_thr, nblocks)
  } else {
    if (length(FDP_thr) == nblocks) {
      FDP_thr_blocks <- FDP_thr
    }
  }
  
  # Computing the total number of iterations
  K <- stability1$params$K + stability2$params$K
  myparams$K <- K
  myparams$seed <- Inf
  
  # Computing selection proportions
  bigstab <- NULL
  if (!is.null(stability1$selprop)) {
    bigstab <- array(NA, dim = dim(stability1$selprop), dimnames = dimnames(stability1$selprop))
    if (graph) {
      for (k in 1:dim(stability1$selprop)[3]) {
        bigstab[, , k] <- (stability1$selprop[, , k] * stability1$params$K + stability2$selprop[, , k] * stability2$params$K) / K
      }
    } else {
      for (k in 1:nrow(stability1$selprop)) {
        bigstab[k, ] <- (stability1$selprop[k, ] * stability1$params$K + stability2$selprop[k, ] * stability2$params$K) / K
      }
    }
  }
  
  # Computing co-membership proportions
  if (inherits(stability1, "clustering")) {
    coprop <- array(NA, dim = dim(stability1$coprop), dimnames = dimnames(stability1$coprop))
    for (k in 1:dim(stability1$coprop)[3]) {
      coprop[, , k] <- (stability1$coprop[, , k] * stability1$params$K + stability2$coprop[, , k] * stability2$params$K) / K
    }
  }
  
  # Concatenating the beta coefficients
  if (inherits(stability1, c("variable_selection", "structural_model", "clustering"))) {
    if (include_beta) {
      Beta <- NULL
      if (!is.null(stability1$Beta)) {
        if (length(dim(stability1$Beta)) == 4) {
          Beta <- array(NA,
                        dim = c(dim(stability1$Beta)[1:2], dim(stability1$Beta)[3] + dim(stability2$Beta)[3], dim(stability1$Beta)[4]),
                        dimnames = list(
                          dimnames(stability1$Beta)[[1]], dimnames(stability1$Beta)[[2]],
                          c(dimnames(stability1$Beta)[[3]], dimnames(stability2$Beta)[[3]]), dimnames(stability1$Beta)[[4]]
                        )
          )
          Beta[, , 1:dim(stability1$Beta)[3], ] <- stability1$Beta
          Beta[, , (dim(stability1$Beta)[3] + 1):dim(Beta)[3], ] <- stability2$Beta
        } else {
          Beta <- array(NA,
                        dim = c(dim(stability1$Beta)[1:2], dim(stability1$Beta)[3] + dim(stability2$Beta)[3]),
                        dimnames = list(
                          dimnames(stability1$Beta)[[1]], dimnames(stability1$Beta)[[2]],
                          c(dimnames(stability1$Beta)[[3]], dimnames(stability2$Beta)[[3]])
                        )
          )
          Beta[, , 1:dim(stability1$Beta)[3]] <- stability1$Beta
          Beta[, , (dim(stability1$Beta)[3] + 1):dim(Beta)[3]] <- stability2$Beta
        }
      }
    }
  }
  
  # Computation of the stability score
  if (inherits(stability1, "graphical_model")) {
    metrics <- StabilityMetrics(
      selprop = bigstab, pk = pk, pi_list = pi_list, K = K, n_cat = n_cat,
      Sequential_template = Sequential_template, graph = graph,
      PFER_method = PFER_method, PFER_thr_blocks = PFER_thr_blocks, FDP_thr_blocks = FDP_thr_blocks
    )
  }
  if (inherits(stability1, c("variable_selection", "structural_model"))) {
    metrics <- StabilityMetrics(
      selprop = bigstab, pk = NULL, pi_list = pi_list, K = K, n_cat = n_cat,
      Sequential_template = Sequential_template, graph = graph,
      PFER_method = PFER_method, PFER_thr_blocks = PFER_thr_blocks, FDP_thr_blocks = FDP_thr_blocks
    )
  }
  if (inherits(stability1, "clustering")) {
    sampled_pairs <- stability1$sampled_pairs + stability2$sampled_pairs
    Sc <- matrix(NA, nrow = dim(coprop)[3], ncol = 1)
    for (k in 1:dim(coprop)[3]) {
      # Clustering on the consensus matrix
      sh_clust <- stats::hclust(stats::as.dist(1 - coprop[, , k]), method = stability1$methods$linkage)
      
      # Identifying stable clusters
      theta <- CoMembership(groups = stats::cutree(sh_clust, k = nc[k]))
      
      Sc[k, 1] <- ConsensusScore(prop = coprop[, , k], K = sampled_pairs, theta = theta)
    }
    Q <- 1 / K * (stability1$params$K * stability1$Q + stability2$params$K * stability2$Q) # weighted average
  }
  
  # Preparing output
  if (inherits(stability1, "graphical_model")) {
    if (nblocks == 1) {
      out <- list(
        S = metrics$S, Lambda = Lambda,
        Q = metrics$Q, Q_s = metrics$Q_s, P = metrics$P,
        PFER = metrics$PFER, FDP = metrics$FDP,
        S_2d = metrics$S_2d, PFER_2d = metrics$PFER_2d, FDP_2d = metrics$FDP_2d,
        selprop = bigstab,
        sign = mysign,
        methods = mymethods,
        params = myparams
      )
    } else {
      out <- list(
        S = metrics$S, Lambda = Lambda,
        Q = metrics$Q, Q_s = metrics$Q_s, P = metrics$P,
        PFER = metrics$PFER, FDP = metrics$FDP,
        S_2d = metrics$S_2d,
        selprop = bigstab,
        sign = mysign,
        methods = mymethods,
        params = myparams
      )
    }
  }
  
  if (inherits(stability1, c("variable_selection", "structural_model"))) {
    if (include_beta) {
      out <- list(
        S = metrics$S, Lambda = Lambda,
        Q = metrics$Q, Q_s = metrics$Q_s, P = metrics$P,
        PFER = metrics$PFER, FDP = metrics$FDP,
        S_2d = metrics$S_2d, PFER_2d = metrics$PFER_2d, FDP_2d = metrics$FDP_2d,
        selprop = bigstab,
        Beta = Beta,
        methods = mymethods,
        params = myparams
      )
    } else {
      out <- list(
        S = metrics$S, Lambda = Lambda,
        Q = metrics$Q, Q_s = metrics$Q_s, P = metrics$P,
        PFER = metrics$PFER, FDP = metrics$FDP,
        S_2d = metrics$S_2d, PFER_2d = metrics$PFER_2d, FDP_2d = metrics$FDP_2d,
        selprop = bigstab,
        methods = mymethods,
        params = myparams
      )
    }
  }
  
  if (inherits(stability1, "clustering")) {
    if (include_beta) {
      out <- list(
        Sc = Sc,
        nc = nc,
        Lambda = Lambda,
        Q = Q,
        coprop = coprop,
        Beta = Beta,
        selprop = bigstab,
        sampled_pairs = sampled_pairs,
        methods = mymethods,
        params = myparams
      )
    } else {
      out <- list(
        Sc = Sc,
        nc = nc,
        Lambda = Lambda,
        Q = Q,
        coprop = coprop,
        selprop = bigstab,
        sampled_pairs = sampled_pairs,
        methods = mymethods,
        params = myparams
      )
    }
  }
  
  # Defining the class
  class(out) <- class(stability1)
  
  return(out)
}

SelectionPerformanceSingle <- function(Asum, cor = NULL, thr = 0.5) {
  # Asum is an adjacency matrix with 3 for TP, 2 for FN, 1 for FP, and 0 for TN
  
  # Preparing objects
  if (is.matrix(Asum)) {
    p <- ncol(Asum)
    N <- p * (p - 1) / 2
    Asum <- Asum[upper.tri(Asum)]
  } else {
    N <- length(Asum)
  }
  
  # Computing the numbers of True/False Positives/Negatives
  TP <- as.numeric(sum(Asum == 3))
  FN <- as.numeric(sum(Asum == 2))
  FP <- as.numeric(sum(Asum == 1))
  TN <- as.numeric(sum(Asum == 0))
  
  # Separation between correlated and independent features based on a threshold in correlation
  if (!is.null(cor)) {
    if (is.matrix(cor)) {
      cor_vect <- cor[upper.tri(cor)]
    } else {
      cor_vect <- cor
    }
    FP_c <- sum((Asum == 1) & (abs(cor_vect) >= thr))
    FP_i <- sum((Asum == 1) & (abs(cor_vect) < thr))
  }
  
  # Computing performances in selection
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  accuracy <- (TP + TN) / N
  if (TP + FP > 0) {
    precision <- TP / (TP + FP)
  } else {
    precision <- 1
  }
  if ((TP + FN) > 0) {
    recall <- TP / (TP + FN)
  } else {
    recall <- 1
  }
  if ((precision > 0) | (recall > 0)) {
    F1_score <- 2 * precision * recall / (precision + recall)
  } else {
    F1_score <- 0
  }
  
  if (is.null(cor)) {
    return(data.frame(
      TP = TP, FN = FN, FP = FP, TN = TN,
      sensitivity = sensitivity, specificity = specificity,
      accuracy = accuracy, precision = precision, recall = recall, F1_score = F1_score
    ))
  } else {
    return(data.frame(
      TP = TP, FN = FN, FP = FP, TN = TN, FP_c = FP_c, FP_i = FP_i,
      sensitivity = sensitivity, specificity = specificity,
      accuracy = accuracy, precision = precision, recall = recall, F1_score = F1_score
    ))
  }
}


SelectionPerformance <- function(theta, theta_star, pk = NULL, cor = NULL, thr = 0.5) {
  # Re-formatting input theta
  if (inherits(theta, c(
    "variable_selection",
    "structural_model",
    "graphical_model",
    "bi_selection"
  ))) {
    if (inherits(theta, c("graphical_model", "structural_model"))) {
      theta <- Adjacency(theta)
    } else {
      if (inherits(theta, "variable_selection")) {
        theta <- SelectedVariables(theta)
        theta <- as.vector(theta)
      } else {
        if ("selected" %in% names(theta)) {
          theta <- theta$selected # PLS
        } else {
          theta <- theta$selectedX # PCA
        }
      }
    }
  }
  
  # Re-formatting input theta_star
  if (inherits(theta_star, c(
    "simulation_regression",
    "simulation_graphical_model",
    "simulation_components",
    "simulation_structural_causal_model"
  ))) {
    theta_star <- theta_star$theta
  }
  if (is.vector(theta)) {
    theta_star <- as.vector(theta_star)
  } else {
    if (ncol(theta) != ncol(theta_star)) {
      theta_star <- theta_star[, 1:ncol(theta)]
    }
  }
  
  # Storing similarities/differences between estimated and true sets
  Asum <- theta + 2 * theta_star
  
  # Extracting block-specific performances
  if (is.null(pk)) {
    if (is.vector(Asum)) {
      out <- SelectionPerformanceSingle(Asum, cor = cor, thr = thr)
    } else {
      if (ncol(theta_star) == nrow(theta_star)) {
        out <- SelectionPerformanceSingle(Asum, cor = cor, thr = thr)
      } else {
        out <- NULL
        for (k in 1:ncol(Asum)) {
          out <- rbind(out, SelectionPerformanceSingle(Asum[, k], cor = cor, thr = thr))
        }
        rownames(out) <- colnames(Asum)
      }
    }
  } else {
    Asum_vect <- Asum[upper.tri(Asum)]
    bigblocks <- BlockMatrix(pk)
    bigblocks_vect <- bigblocks[upper.tri(bigblocks)]
    if (!is.null(cor)) {
      cor_vect <- cor[upper.tri(cor)]
    } else {
      cor_vect <- NULL
    }
    
    out <- SelectionPerformanceSingle(Asum, cor = cor, thr = thr)
    for (k in sort(unique(bigblocks_vect))) {
      tmp <- SelectionPerformanceSingle(Asum_vect[bigblocks_vect == k],
                                        cor = cor_vect[bigblocks_vect == k], thr = thr
      )
      out <- rbind(out, tmp)
    }
  }
  
  return(out)
}


Stable <- function(stability, argmax_id = NULL, linkage = "complete") {
  if (inherits(stability, c("variable_selection", "bi_selection"))) {
    out <- SelectedVariables(stability = stability, argmax_id = argmax_id)
  }
  
  if (inherits(stability, c("graphical_model", "structural_model"))) {
    out <- Adjacency(stability = stability, argmax_id = argmax_id)
  }
  
  if (inherits(stability, c("clustering"))) {
    out <- Clusters(stability = stability, linkage = linkage, argmax_id = argmax_id)
  }
  
  return(out)
}

SelectedVariables <- function(stability, argmax_id = NULL) {
  if (!inherits(stability, c("variable_selection", "structural_model", "bi_selection"))) {
    stop("Invalid input for argument 'stability'. This function only applies to outputs from VariableSelection() or BiSelection().")
  }
  
  if (inherits(stability, c("variable_selection", "structural_model"))) {
    if (is.null(argmax_id)) {
      argmax_id <- ArgmaxId(stability)
      argmax <- Argmax(stability)
    } else {
      argmax <- c(NA, stability$params$pi_list[argmax_id[2]])
    }
    stability_selected <- ifelse(stability$selprop[argmax_id[1], ] >= argmax[2],
                                 yes = 1, no = 0
    )
  }
  
  if (inherits(stability, "bi_selection")) {
    if (is.null(argmax_id)) {
      stability_selected <- stability$selectedX
    } else {
      stop("Invalid input for argument 'argmax_id'. Arbitrary choice of parameters is not supported for the output of BiSelection().")
    }
  }
  
  return(stability_selected)
}

SelectionProportionsRegression <- function(stability, argmax_id = NULL) {
  if (is.null(argmax_id)) {
    argmax_id <- ArgmaxId(stability)
  }
  m <- stability$selprop[argmax_id[1], ]
  return(m)
}

SelectionProportions <- function(stability, argmax_id = NULL) {
  out <- NULL
  
  if (inherits(stability, "graphical_model")) {
    out <- SelectionProportionsGraphical(stability = stability, argmax_id = argmax_id)
  }
  if (inherits(stability, "variable_selection")) {
    out <- SelectionProportionsRegression(stability = stability, argmax_id = argmax_id)
  }
  if (inherits(stability, "clustering")) {
    argmax_id <- ArgmaxId(stability)
    out <- SelectionProportionsRegression(stability = stability, argmax_id = argmax_id)
  }
  if (inherits(stability, "bi_selection")) {
    out <- stability$selpropX
  }
  
  return(out)
}

ArgmaxId <- function(stability = NULL, S = NULL) {
  if ((is.null(stability)) & (is.null(S))) {
    stop("Invalid input. One of the two arguments has to be specified: 'stability' or 'S'.")
  }
  
  clustering <- ifelse(inherits(stability, "clustering"), yes = TRUE, no = FALSE)
  
  if (is.null(S)) {
    if (clustering) {
      # If multiple solutions, prioritising many clusters over small number of features
      Sc <- round(stability$Sc, digits = 4)
      argmax_id <- which(Sc == max(Sc, na.rm = TRUE))
      argmax_id <- argmax_id[which(stability$nc[argmax_id] == max(stability$nc[argmax_id]))]
      argmax_id <- cbind(min(argmax_id))
    } else {
      argmax_id <- matrix(NA, nrow = ncol(stability$Lambda), ncol = 2)
      if (is.null(stability$params$lambda_other_blocks) & (length(stability$params$pk) > 1)) {
        id <- which.max(apply(stability$S, 1, sum, na.rm = TRUE))
        argmax_id[, 1] <- rep(id, nrow(argmax_id))
        for (block_id in 1:ncol(stability$Lambda)) {
          if (!is.na(stability$P[id, block_id])) {
            argmax_id[block_id, 2] <- which(stability$params$pi_list == stability$P[id, block_id])
          }
        }
      } else {
        for (block_id in 1:ncol(stability$Lambda)) {
          if (ncol(stability$Lambda) == 1) {
            myS <- stability$S
          } else {
            myS <- stability$S[, block_id, drop = FALSE]
          }
          myS[is.na(myS)] <- 0
          myid <- which.max(myS[, 1])
          argmax_id[block_id, ] <- c(myid, which(stability$params$pi_list == stability$P[myid, block_id]))
        }
      }
    }
  } else {
    argmax_id <- matrix(NA, nrow = 1, ncol = 2)
    myS <- apply(S, 1, max, na.rm = TRUE)
    myS[is.na(myS)] <- 0
    myid <- which.max(myS)
    argmax_id[1, ] <- c(myid, max(which(S[myid, ] == myS[myid])))
  }
  if (clustering) {
    colnames(argmax_id) <- c("row_id")
  } else {
    colnames(argmax_id) <- c("lambda_id", "pi_id")
  }
  return(argmax_id)
}

Argmax <- function(stability) {
  clustering <- ifelse(inherits(stability, "clustering"), yes = TRUE, no = FALSE)
  if (inherits(stability, "bi_selection")) {
    argmax <- stability$summary
    argmax <- argmax[, colnames(argmax) != "S", drop = FALSE]
  } else {
    if (clustering) {
      id <- ArgmaxId(stability = stability)
      argmax <- matrix(
        c(
          stability$nc[id[1], 1],
          stability$Lambda[id[1], 1]
        ),
        ncol = 2
      )
    } else {
      argmax <- matrix(NA, nrow = ncol(stability$Lambda), ncol = 2)
      if (is.null(stability$params$lambda_other_blocks) & (length(stability$params$pk) > 1)) {
        id <- which.max(apply(stability$S, 1, sum, na.rm = TRUE))
        argmax[, 1] <- stability$Lambda[id, ]
        argmax[, 2] <- stability$P[id, ]
      } else {
        for (block_id in 1:ncol(stability$Lambda)) {
          if (ncol(stability$Lambda) == 1) {
            myS <- stability$S
          } else {
            myS <- stability$S[, block_id, drop = FALSE]
          }
          myS[is.na(myS)] <- 0
          myid <- which.max(myS[, 1])
          argmax[block_id, ] <- c(stability$Lambda[myid, block_id], stability$P[myid, block_id])
        }
      }
    }
    if (clustering) {
      colnames(argmax) <- c("nc", "lambda")
    } else {
      colnames(argmax) <- c("lambda", "pi")
    }
  }
  
  return(argmax)
}

PFER <- function(q, pi, N, K, PFER_method = "MB") {
  # Checking the inputs (PFER_method)
  PFER_method <- as.character(PFER_method)
  if ((length(PFER_method) != 1) | (!PFER_method %in% c("MB", "SS"))) {
    stop("Invalid input for argument 'PFER_method'. Possible values are: 'MB' or 'SS'.")
  }
  
  if (pi > 0.5) {
    # Computing upper-bound of the PFER using approach proposed by MB
    if (PFER_method == "MB") {
      upperbound <- 1 / (2 * pi - 1) * q^2 / N
    }
    
    # Computing upper-bound of the PFER using approach proposed by SS
    if (PFER_method == "SS") {
      cutoff <- pi
      B <- ceiling(K / 2)
      theta <- q / N
      if (cutoff <= 3 / 4) {
        tmp <- 2 * (2 * cutoff - 1 - 1 / (2 * B))
      } else {
        tmp <- (1 + 1 / B) / (4 * (1 - cutoff + 1 / (2 * B)))
      }
      upperbound <- q^2 / N / tmp
      
      # Setting to Inf if "out of bounds"
      if ((cutoff < 1 / 2 + min(theta^2, 1 / (2 * B) + 3 / 4 * theta^2)) | (cutoff > 1)) {
        upperbound <- Inf
      }
    }
  } else {
    upperbound <- Inf
  }
  
  # Re-formatting the upperbound
  if (is.na(upperbound)) {
    upperbound <- Inf
  }
  
  return(upperbound)
}

FDP <- function(selprop, PFER, pi) {
  # Preparing objects
  if (is.matrix(selprop)) {
    selprop <- selprop[upper.tri(selprop)]
  }
  
  # Computing the number of stable edges
  S <- sum(selprop >= pi, na.rm = TRUE)
  
  # Computing the proportion of false discoveries among discoveries (False Discovery Proportion)
  if (S != 0) {
    FDP <- PFER / S
  } else {
    FDP <- 0
  }
  
  return(FDP)
}

LambdaSequence <- function(lmax, lmin, cardinal = 100) {
  return(exp(seq(log(lmax), log(lmin), length.out = cardinal)))
  # return(seq(sqrt(lmax),sqrt(lmin),length.out=cardinal)^2)
}

Resample <- function(data, family = NULL, tau = 0.5, resampling = "subsampling", ...) {
  # Preparing the data
  if (is.factor(data)) {
    data <- as.character(factor(data, levels = levels(data), labels = seq(1, length(levels(data))) - 1))
  }
  if (is.vector(data)) {
    data <- matrix(data, ncol = 1)
  }
  if (!is.null(family)) {
    if (family == "multinomial") {
      if (is.matrix(data)) {
        data <- DummyToCategories(x = data, verbose = FALSE)
      }
    }
  }
  
  # if (!resampling %in% c("subsampling", "bootstrap")) {
  if (is.function(resampling)) {
    # s <- do.call(get(resampling), args = list(data = data, tau = tau, ...))
    s <- do.call(resampling, args = list(data = data, tau = tau, ...))
  } else {
    if (!resampling %in% c("subsampling", "bootstrap")) {
      stop("Invalid input for argument 'resampling'. It must be a function or a character string: 'subsampling' or 'bootstrap'.")
    } else {
      # Using or not replacement in resampling
      replacement <- ifelse(resampling == "subsampling", yes = FALSE, no = TRUE)
      
      # Definition of the size of sub/bootstrap sample
      if (replacement) {
        tau <- 1
      }
      
      # Resampling procedure
      if (!is.null(family)) {
        # Resampling for regression models
        if (family %in% c("gaussian", "poisson", "mgaussian")) {
          s <- sample(nrow(data), size = tau * nrow(data), replace = replacement)
        }
        if (family == "binomial") {
          if (ncol(data) > 1) {
            data <- cbind(apply(data, 1, sum)) # to ensure balanced classes for PLS-DA
          }
          s <- NULL
          for (mycat in levels(factor(data))) {
            scat <- sample(which(data == mycat), size = tau * sum(data == mycat), replace = replacement)
            s <- c(s, scat)
          }
        }
        if (family == "multinomial") {
          s <- NULL
          for (mycat in levels(factor(data))) {
            scat <- sample(which(data == mycat), size = tau * sum(data == mycat), replace = replacement)
            s <- c(s, scat)
          }
        }
        if (family == "cox") {
          s0 <- sample(which(data[, 2] == "0"), size = tau * sum(data[, 2] == "0"), replace = replacement)
          s1 <- sample(which(data[, 2] == "1"), size = tau * sum(data[, 2] == "1"), replace = replacement)
          s <- c(s0, s1)
        }
      } else {
        # Resampling for network models
        s <- sample(1:nrow(data), size = tau * nrow(data), replace = replacement)
      }
    }
  }
  return(s)
}

ConsensusScore <- function(prop, K, theta) {
  # Calculating the within and between sums
  X_w <- sum(prop * K * theta) 
  X_b <- sum(prop * K * (1 - theta))
  N_w <- sum(K * theta)
  N_b <- sum(K * (1 - theta))
  
  # Calculating the z statistic 
  p_w <- X_w / N_w 
  p_b <- X_b / N_b 
  p_0 <- (X_w + X_b) / (N_w + N_b)
  score <- (p_w - p_b) / sqrt(p_0 * (1 - p_0) * (1 / N_w + 1 / N_b))

  # # Calculating consensus score as negative log p-value
  # score <- -stats::pnorm(z, lower.tail = FALSE, log.p = TRUE)
  
  return(score)
}

BinomialProbabilities <- function(q, N, pi, K, n_cat = 3) {
  if (n_cat == 2) {
    # Definition of the threshold in selection counts
    thr <- round(K * pi) # Threshold above (>=) which the feature is stably selected
    
    # Probability of observing a selection count below thr_down under the null (uniform selection)
    p_0 <- stats::pbinom(thr - 1, size = K, prob = q / N, log.p = TRUE) # proportion < pi
    
    # Probability of observing a selection count above thr_up under the null
    p_1 <- stats::pbinom(thr - 1, size = K, prob = q / N, lower.tail = FALSE, log.p = TRUE) # proportion >= pi
    
    # Checking consistency between the three computed probabilities (should sum to 1)
    if (abs(exp(p_0) + exp(p_1) - 1) > 1e-3) {
      message(paste("N:", N))
      message(paste("q:", q))
      message(paste("K:", K))
      message(paste("pi:", pi))
      stop(paste0("Probabilities do not sum to 1 (Binomial distribution) \n p_0+p_1=", exp(p_0) + exp(p_1)))
    }
    
    # Output the two probabilities under the assumption of uniform selection procedure
    return(list(p_0 = p_0, p_1 = p_1))
  }
  
  if (n_cat == 3) {
    # Definition of the two thresholds in selection counts
    thr_down <- round(K * (1 - pi)) # Threshold below (<=) which the feature is stable-out
    thr_up <- round(K * pi) # Threshold above (>=) which the feature is stable-in
    
    # Probability of observing a selection count below thr_down under the null (uniform selection)
    p_1 <- stats::pbinom(thr_down, size = K, prob = q / N, log.p = TRUE) # proportion <= (1-pi)
    
    # Probability of observing a selection count between thr_down and thr_up under the null
    if ((thr_down) >= (thr_up - 1)) {
      # Not possible to compute (i.e. negative number of unstable features)
      p_2 <- NA
    } else {
      # Using cumulative probabilities
      p_2 <- log(stats::pbinom(thr_up - 1, size = K, prob = q / N) - stats::pbinom(thr_down, size = K, prob = q / N)) # 1-pi < proportion < pi
      
      # Using sum of probabilities (should not be necessary)
      if (is.infinite(p_2) | is.na(p_2)) {
        p_2 <- 0
        for (i in seq(thr_down + 1, thr_up - 1)) {
          p_2 <- p_2 + stats::dbinom(i, size = K, prob = q / N)
        }
        p_2 <- log(p_2)
      }
    }
    
    # Probability of observing a selection count above thr_up under the null
    p_3 <- stats::pbinom(thr_up - 1, size = K, prob = q / N, lower.tail = FALSE, log.p = TRUE) # proportion >= pi
    
    # Checking consistency between the three computed probabilities (should sum to 1)
    if (!is.na(p_2)) {
      if (abs(exp(p_1) + exp(p_2) + exp(p_3) - 1) > 1e-3) {
        message(paste("N:", N))
        message(paste("q:", q))
        message(paste("K:", K))
        message(paste("pi:", pi))
        stop(paste0("Probabilities do not sum to 1 (Binomial distribution) \n p_1+p_2+p_3=", exp(p_1) + exp(p_2) + exp(p_3)))
      }
    }
    
    # Output the three probabilities under the assumption of uniform selection procedure
    return(list(p_1 = p_1, p_2 = p_2, p_3 = p_3))
  }
}

StabilityScore <- function(selprop, pi_list = seq(0.6, 0.9, by = 0.01), K, n_cat = 3, group = NULL) {
  # Preparing objects
  if (is.matrix(selprop)) {
    selprop <- selprop[upper.tri(selprop)]
  }
  
  # Using group penalisation (extracting one per group)
  if (!is.null(group)) {
    selprop <- selprop[cumsum(group)]
  }
  
  # Computing the number of features (edges/variables)
  N <- sum(!is.na(selprop))
  
  # Computing the average number of selected features
  q <- round(sum(selprop, na.rm = TRUE))
  
  # Loop over the values of pi
  score <- rep(NA, length(pi_list))
  for (i in 1:length(pi_list)) {
    pi <- pi_list[i]
    
    # Computing the probabilities of being stable-in, stable-out or unstable under the null (uniform selection) 
    p_vect <- BinomialProbabilities(q, N, pi, K, n_cat = n_cat) 
    
    # Computing the log-likelihood
    if (any(is.na(p_vect))) {
      # Returning NA if not possible to compute (e.g. negative number of unstable features, as with pi<=0.5)
      l <- NA
    } else {
      if (n_cat == 2) {
        S_0 <- sum(selprop < pi, na.rm = TRUE) # Number of not stably selected features
        S_1 <- sum(selprop >= pi, na.rm = TRUE) # Number of stably selected features
        
        # Checking consistency
        if (S_0 + S_1 != N) {
          stop(paste0("Inconsistency in number of edges \n S_0+S_1=", S_0 + S_1, " instead of ", N))
        }
        
        # Log-likelihood
        l <- S_0 * p_vect$p_0 + S_1 * p_vect$p_1
      }
      
      if (n_cat == 3) {
        S_0 <- sum(selprop <= (1 - pi), na.rm = TRUE) # Number of stable-out features
        S_1 <- sum(selprop >= pi, na.rm = TRUE) # Number of stable-in features
        U <- sum((selprop < pi) & (selprop > (1 - pi)), na.rm = TRUE) # Number of unstable features
        
        # Checking consistency
        if (S_0 + S_1 + U != N) {
          stop(paste0("Inconsistency in number of edges \n S_0+S_1+U=", S_0 + S_1 + U, " instead of ", N))
        }
        
        # Log-likelihood
        l <- S_0 * p_vect$p_1 + U * p_vect$p_2 + S_1 * p_vect$p_3
      }
      
      # Re-formatting if infinite
      if (is.infinite(l)) {
        l <- NA
      }
    }
    
    # Getting the stability score
    score[i] <- -l
  }
  
  return(score)
}

StabilityMetrics <- function(selprop, pk = NULL, pi_list = seq(0.6, 0.9, by = 0.01),
                             K = 100, n_cat = NULL,
                             PFER_method = "MB", PFER_thr_blocks = Inf, FDP_thr_blocks = Inf,
                             Sequential_template = NULL, graph = TRUE, group = NULL) {
  if (graph) {
    nlambda <- dim(selprop)[3]
  } else {
    nlambda <- nrow(selprop)
  }
  
  # Extracting pk
  if (is.null(pk)) {
    pk <- ncol(selprop)
  }
  
  if (is.null(Sequential_template)) {
    Sequential_template <- matrix(TRUE, nrow = nlambda, ncol = 1)
  }
  
  # Creating matrix with block indices
  nblocks <- 1
  if (graph) { # to avoid memory issues in high dimensional variable selection
    bigblocks <- BlockMatrix(pk)
    nblocks <- length(pk) * (length(pk) + 1) / 2
    bigblocks_vect <- factor(bigblocks[upper.tri(bigblocks)], levels = 1:nblocks)
    N_blocks <- unname(table(bigblocks_vect))
    blocks <- levels(bigblocks_vect)
    names(N_blocks) <- blocks
  }
  
  # Initialising objects to be filled
  Q <- Q_s <- P <- matrix(NA, nrow = nlambda, ncol = nblocks)
  best_loglik <- best_PFER <- best_FDP <- matrix(NA, nrow = nlambda, ncol = nblocks)
  if (nblocks == 1) {
    loglik <- PFER <- FDP <- matrix(NA, ncol = length(pi_list), nrow = nlambda)
  } else {
    loglik <- array(NA, dim = c(nlambda, length(pi_list), nblocks))
  }
  
  # Computing the metrics for each value of lambda
  for (k in 1:nlambda) {
    # Extracting corresponding selection proportions
    if (graph) {
      stab_iter <- selprop[, , k]
    } else {
      stab_iter <- selprop[k, ]
    }
    
    # Computing stability score with block-specific pi
    for (block_id in 1:nblocks) {
      if (Sequential_template[k, block_id]) {
        if (graph) {
          stab_iter_block <- stab_iter[(bigblocks == block_id) & (upper.tri(bigblocks))] # selection proportions in the block
        } else {
          stab_iter_block <- stab_iter
        }
        
        # Using group penalisation (extracting one per group)
        if (!is.null(group)) {
          stab_iter_block <- stab_iter_block[cumsum(group)]
        }
        
        q_block <- round(sum(stab_iter_block, na.rm = TRUE)) # average number of edges selected by the original procedure in the block
        Q[k, block_id] <- q_block
        N_block <- length(stab_iter_block) # maximum number of edges in the block
        tmp_loglik <- tmp_PFERs <- tmp_FDPs <- rep(NA, length(pi_list))
        
        # Computing error rates and stability score for different values of pi
        for (j in 1:length(pi_list)) {
          pi <- pi_list[j]
          tmp_PFERs[j] <- PFER(q = q_block, pi = pi, N = N_block, K = K, PFER_method = PFER_method)
          tmp_FDPs[j] <- FDP(selprop = stab_iter_block, PFER = tmp_PFERs[j], pi = pi)
          if ((tmp_PFERs[j] <= PFER_thr_blocks[block_id]) & (tmp_FDPs[j] <= FDP_thr_blocks[block_id])) {
            # Computing stability score (group penalisation is accounted for above so no need here)
            if (is.null(n_cat)) {
              theta <- ifelse(stab_iter_block >= pi, yes = 1, no = 0)
              tmp_loglik[j] <- ConsensusScore(prop = stab_iter_block, K = K, theta = theta)
            } else {
              tmp_loglik[j] <- StabilityScore(selprop = stab_iter_block, pi_list = pi, K = K, n_cat = n_cat, group = NULL)
            }
          }
        }
        
        # Storing stability score in a matrix if only one block
        if (nblocks == 1) {
          loglik[k, ] <- tmp_loglik
          PFER[k, ] <- tmp_PFERs
          FDP[k, ] <- tmp_FDPs
        } else {
          loglik[k, , block_id] <- tmp_loglik
        }
        
        # Keeping best stability score and other parameters at the max
        if (any(!is.na(tmp_loglik))) {
          tmp_loglik[is.na(tmp_loglik)] <- 0
          
          # Extracting the best pi (highest threshold if multiple choices)
          myid <- which(tmp_loglik == max(tmp_loglik, na.rm = TRUE))
          myid <- max(myid)
          
          # Extracting corresponding metrics
          tmp_loglik[which(tmp_loglik == 0)] <- NA
          best_loglik[k, block_id] <- tmp_loglik[myid]
          P[k, block_id] <- pi_list[myid]
          Q_s[k, block_id] <- sum(stab_iter_block >= pi_list[myid], na.rm = TRUE)
          best_PFER[k, block_id] <- tmp_PFERs[myid]
          best_FDP[k, block_id] <- tmp_FDPs[myid]
        }
      }
    }
  }
  best_loglik_blocks <- best_loglik
  best_loglik <- matrix(apply(best_loglik, 1, sum), ncol = 1)
  
  if (nblocks == 1) {
    return(list(
      S = best_loglik_blocks,
      Q = Q, Q_s = Q_s, P = P,
      PFER = best_PFER, FDP = best_FDP,
      S_2d = loglik, PFER_2d = PFER, FDP_2d = FDP
    ))
  } else {
    return(list(
      S = best_loglik_blocks,
      Q = Q, Q_s = Q_s, P = P,
      PFER = best_PFER, FDP = best_FDP,
      S_2d = loglik
    ))
  }
}

PenalisedRegression <- function(xdata, ydata, Lambda = NULL, family,
                                penalisation = c("classic", "randomised", "adaptive"),
                                gamma = NULL,
                                # model = c("glmnet","glinternet"), 
                                ...) {
  # Checking that input data are matrices
  xdata <- as.matrix(xdata)
  ydata <- as.matrix(ydata)
  
  # Defining the type of penalised regression
  penalisation <- match.arg(penalisation)
  
  # Setting the default gamma
  if (is.null(gamma)) {
    gamma <- switch(penalisation,
                    classic = NULL,
                    randomised = 0.5,
                    adaptive = 2
    )
  }
  
  # Storing extra arguments
  extra_args <- list(...)
  
  # Preparing the data for randomised lasso
  if (penalisation == "randomised") {
    xdata <- t(t(xdata) * stats::runif(ncol(xdata), min = gamma, max = 1))
  }
  
  # Preparing the data for adaptive lasso
  if (penalisation == "adaptive") {
    # Running unpenalised model to get the weights
    if (family == "multinomial") {
      # Extracting relevant extra arguments (excluding penalty.factor here)
      tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = glmnet::glmnet)
      tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("x", "y", "lambda", "standardize", "family", "type.multinomial", "penalty.factor")]
      
      # Running model
      myweights <- stats::coef(do.call(glmnet::glmnet, args = c(
        list(
          x = xdata,
          y = ydata,
          family = family,
          lambda = 0,
          standardize = FALSE,
          # thresh = 1e-15,
          type.multinomial = "grouped"
        ),
        tmp_extra_args
      )))[-1, 1]
    } else {
      # Extracting relevant extra arguments (excluding penalty.factor here)
      tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = glmnet::glmnet)
      tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("x", "y", "lambda", "family", "penalty.factor")]
      
      # Running model
      myweights <- stats::coef(do.call(glmnet::glmnet, args = c(
        list(
          x = xdata,
          y = ydata,
          lambda = 0,
          standardize = FALSE,
          # thresh = 1e-15, # could cause convergence issues
          family = family
        ),
        tmp_extra_args
      )))[-1, 1]
    }
    
    # Using the weights as penalty factors
    myweights <- 1 / (abs(myweights)^gamma)
    if ("penalty.factor" %in% names(extra_args)) {
      extra_args$penalty.factor <- extra_args$penalty.factor * myweights
    } else {
      extra_args$penalty.factor <- myweights
    }
  }
  
  # Running the regression
  if (family == "multinomial") {
    # Extracting relevant extra arguments
    tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = glmnet::glmnet)
    tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("x", "y", "lambda", "standardize", "family", "type.multinomial")]
    
    # Running model
    mymodel <- do.call(glmnet::glmnet, args = c(
      list(
        x = xdata,
        y = ydata,
        family = family,
        lambda = Lambda,
        standardize = FALSE,
        type.multinomial = "grouped"
      ),
      tmp_extra_args
    ))
  } else {
    # Extracting relevant extra arguments
    tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = glmnet::glmnet)
    tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("x", "y", "lambda", "standardize", "family")]
    
    # Running model
    mymodel <- do.call(glmnet::glmnet, args = c(
      list(
        x = xdata,
        y = ydata,
        lambda = Lambda,
        standardize = FALSE,
        family = family
      ),
      tmp_extra_args
    ))
  }
  if (!is.infinite(mymodel$lambda[1])) {
    # Extracting and formatting the beta coefficients
    if (!family %in% c("mgaussian", "multinomial")) {
      if (length(Lambda) > 1) {
        mybeta <- suppressWarnings({
          do.call(stats::coef,
                  args = c(list(object = mymodel, s = Lambda, exact = TRUE, x = xdata, y = ydata), tmp_extra_args)
          )
        })
      } else {
        mybeta <- suppressWarnings({
          stats::coef(mymodel, s = Lambda)
        })
      }
      mybeta <- t(as.matrix(mybeta))
      
      # Preparing the outputs
      beta_full <- mybeta[, colnames(xdata), drop = FALSE] # removing the intercept if included
      if ("penalty.factor" %in% names(extra_args)) {
        selected <- ifelse(mybeta[, colnames(xdata)[which(extra_args$penalty.factor != 0)], drop = FALSE] != 0, yes = 1, no = 0)
      } else {
        selected <- ifelse(mybeta[, colnames(xdata), drop = FALSE] != 0, yes = 1, no = 0)
      }
    } else {
      if (family == "mgaussian") {
        mybeta <- array(NA,
                        dim = c(length(Lambda), ncol(xdata), ncol(ydata)),
                        dimnames = list(1:length(Lambda), colnames(xdata), colnames(ydata))
        )
        if (length(Lambda) > 1) {
          tmpcoefs <- suppressWarnings({
            do.call(stats::coef,
                    args = c(list(object = mymodel, s = Lambda, exact = TRUE, x = xdata, y = ydata), tmp_extra_args)
            )
          })
        } else {
          tmpcoefs <- suppressWarnings({
            stats::coef(mymodel, s = Lambda)
          })
        }
        for (y_id in 1:ncol(ydata)) {
          tmpbeta <- tmpcoefs[[y_id]]
          tmpbeta <- t(as.matrix(tmpbeta))
          tmpbeta <- tmpbeta[, colnames(xdata), drop = FALSE] # removing the intercept if included
          mybeta[rownames(tmpbeta), colnames(tmpbeta), y_id] <- tmpbeta
        }
      }
      if (family == "multinomial") {
        y_levels <- sort(unique(ydata))
        mybeta <- array(NA,
                        dim = c(length(Lambda), ncol(xdata), length(y_levels)),
                        dimnames = list(1:length(Lambda), colnames(xdata), paste0("Y", y_levels))
        )
        if (length(Lambda) > 1) {
          tmpcoefs <- suppressWarnings({
            do.call(stats::coef,
                    args = c(list(object = mymodel, s = Lambda, exact = TRUE, x = xdata, y = ydata), tmp_extra_args)
            )
          })
        } else {
          tmpcoefs <- suppressWarnings({
            stats::coef(mymodel, s = Lambda)
          })
        }
        for (y_id in 1:length(y_levels)) {
          tmpbeta <- tmpcoefs[[y_id]]
          tmpbeta <- t(as.matrix(tmpbeta))
          tmpbeta <- tmpbeta[, colnames(xdata), drop = FALSE] # removing the intercept if included
          mybeta[rownames(tmpbeta), colnames(tmpbeta), y_id] <- tmpbeta
        }
      }
      
      # Preparing the outputs
      if ("penalty.factor" %in% names(extra_args)) {
        selected <- ifelse(abind::adrop(mybeta[, colnames(xdata)[which(extra_args$penalty.factor != 0)], 1, drop = FALSE], drop = 3) != 0, yes = 1, no = 0)
      } else {
        selected <- ifelse(abind::adrop(mybeta[, , 1, drop = FALSE], drop = 3) != 0, yes = 1, no = 0)
      }
      beta_full <- mybeta
    }
  } else {
    # Returning infinite beta is the model failed
    selected <- beta_full <- Inf
  }
  
  return(list(selected = selected, beta_full = beta_full))
}

SelectionAlgo <- function(xdata, ydata = NULL,
                          Lambda, group_x = NULL, scale = TRUE,
                          family = NULL,
                          implementation = PenalisedRegression, ...) {
  # Making sure none of the variables has a null standard deviation
  mysd <- rep(NA, ncol(xdata))
  for (j in 1:ncol(xdata)) {
    mysd[j] <- stats::sd(xdata[, j])
  }
  names(mysd) <- colnames(xdata)
  if (any(mysd == 0)) {
    for (k in which(mysd == 0)) {
      xdata[, k] <- xdata[, k] + stats::rnorm(n = nrow(xdata), sd = min(mysd[mysd != 0]) / 100)
    }
  }
  
  # Scaling the predictor data
  if (scale) {
    xdata <- scale(xdata)
  }

  # Applying user-defined function for variable selection
  mybeta <- do.call(implementation, args = list(xdata = xdata, ydata = ydata, Lambda = Lambda, group_x = group_x, family = family))
  selected <- mybeta$selected
  beta_full <- mybeta$beta_full
  
  # Checking row and column names
  if (is.null(rownames(selected)) | is.null(rownames(beta_full))) {
    rownames(selected) <- rownames(beta_full) <- paste0("s", seq(0, nrow(beta_full) - 1))
  }
  if (is.null(colnames(selected)) | is.null(colnames(beta_full))) {
    colnames(selected) <- colnames(beta_full) <- paste0("coef", seq(1, ncol(beta_full)))
  }
  
  # Setting the beta coefficient to zero for predictors with always the same value (null standard deviation)
  if (!is.infinite(selected[1])) {
    if (any(mysd == 0)) {
      ids_no_sd <- intersect(names(mysd)[which(mysd == 0)], colnames(selected))
      if (length(ids_no_sd) > 0) {
        selected[, ids_no_sd] <- 0
        if (length(dim(beta_full)) == 2) {
          beta_full[, ids_no_sd] <- 0
        }
        if (length(dim(beta_full)) == 3) {
          beta_full[, ids_no_sd, ] <- 0
        }
      }
    }
  }
  
  return(list(selected = selected, beta_full = beta_full))
}
SerialRegression <- function(xdata, ydata = NULL, Lambda, pi_list = seq(0.6, 0.9, by = 0.01),
                             K = 100, tau = 0.5, seed = 1, n_cat = 3,
                             family = "gaussian", implementation = PenalisedRegression,
                             resampling = "subsampling", cpss = FALSE,
                             PFER_method = "MB", PFER_thr = Inf, FDP_thr = Inf,
                             group_x = NULL, group_penalisation = FALSE,
                             output_data = FALSE, verbose = TRUE, ...) {
  # Using complementary pairs for SS
  if (PFER_method == "SS") {
    cpss <- TRUE
  }
  
  # Defining K if using complementary pairs
  if (cpss) {
    K <- ceiling(K / 2) * 2
    tau <- 0.5
  }
  
  # Initialising objects to be filled
  N <- N_block <- ncol(xdata)
  
  # Initialising the arrays
  s <- Resample(data = ydata, family = family, tau = tau, resampling = resampling, ...)
  Xsub <- xdata[s, , drop = FALSE]
  Ysub <- ydata[s, , drop = FALSE]
  mybeta <- SelectionAlgo(
    xdata = Xsub, ydata = Ysub,
    Lambda = Lambda[, 1], group_x = group_x,
    family = family, implementation = implementation, ...
  )
  Beta <- array(0, dim = c(nrow(mybeta$selected), ncol(mybeta$selected), K))
  rownames(Beta) <- rownames(mybeta$selected)
  colnames(Beta) <- colnames(mybeta$selected)
  if (length(dim(mybeta$beta_full)) == 2) {
    Beta_full <- array(0,
                       dim = c(dim(mybeta$beta_full)[1], dim(mybeta$beta_full)[2], K),
                       dimnames = list(rownames(mybeta$beta_full), dimnames(mybeta$beta_full)[[2]], NULL)
    )
  } else {
    if (length(dim(mybeta$beta_full)) == 3) {
      Beta_full <- array(0,
                         dim = c(dim(mybeta$beta_full)[1], dim(mybeta$beta_full)[2], K, dim(mybeta$beta_full)[3]),
                         dimnames = list(rownames(mybeta$beta_full), dimnames(mybeta$beta_full)[[2]], NULL, dimnames(mybeta$beta_full)[[3]])
      )
    } else {
      stop(paste0("Invalid output from the variable selection function provided in 'implementation'. The output 'beta_full' must be an array with 2 or 3 dimensions."))
    }
  }
  
  # Setting seed for reproducibility
  withr::local_seed(seed)
  
  # Computation of the selection proportions over Lambda
  if (verbose) {
    pb <- utils::txtProgressBar(style = 3)
  }
  if (!cpss) {
    for (k in 1:K) {
      s <- Resample(data = ydata, family = family, tau = tau, resampling = resampling, ...)
      Xsub <- xdata[s, , drop = FALSE]
      Ysub <- ydata[s, , drop = FALSE]
      mybeta <- SelectionAlgo(
        xdata = Xsub, ydata = Ysub,
        Lambda = Lambda[, 1], group_x = group_x,
        family = family, implementation = implementation, ...
      )
      
      # Resampling if model failed to converge
      while (is.infinite(mybeta$selected[1])) {
        s <- Resample(data = ydata, family = family, tau = tau, resampling = resampling, ...)
        Xsub <- xdata[s, , drop = FALSE]
        Ysub <- ydata[s, , drop = FALSE]
        mybeta <- SelectionAlgo(
          xdata = Xsub, ydata = Ysub,
          Lambda = Lambda[, 1], group_x = group_x,
          family = family, implementation = implementation, ...
        )
      }
      
      # Storing (one set of) beta coefficients, used to define set of selected variables
      Beta[rownames(mybeta$selected), colnames(mybeta$selected), k] <- mybeta$selected
      
      # Storing all beta coefficients
      if (length(dim(Beta_full)) == 3) {
        Beta_full[rownames(mybeta$beta_full), colnames(mybeta$beta_full), k] <- mybeta$beta_full
      } else {
        Beta_full[rownames(mybeta$beta_full), colnames(mybeta$beta_full), k, ] <- mybeta$beta_full
      }
      
      if (verbose) {
        utils::setTxtProgressBar(pb, k / K)
      }
    }
    
    # Computing the selection proportions
    bigstab <- matrix(NA, nrow = nrow(Beta), ncol = ncol(Beta))
    colnames(bigstab) <- colnames(Beta)
    rownames(bigstab) <- rownames(Beta)
    for (i in 1:nrow(Beta)) {
      for (j in 1:ncol(Beta)) {
        bigstab[i, j] <- sum(Beta[i, j, ] != 0, na.rm = TRUE) / sum(!is.na(Beta[i, j, ]))
      }
    }
  } else {
    for (k in 1:ceiling(K / 2)) {
      s <- Resample(data = ydata, family = family, tau = tau, resampling = resampling, ...)
      
      # First subset
      Xsub <- xdata[s, , drop = FALSE]
      Ysub <- ydata[s, , drop = FALSE]
      mybeta1 <- SelectionAlgo(
        xdata = Xsub, ydata = Ysub,
        Lambda = Lambda[, 1], group_x = group_x,
        family = family, implementation = implementation, ...
      )
      
      # Complementary subset
      Xsub <- xdata[seq(1, nrow(xdata))[!seq(1, nrow(xdata)) %in% s], , drop = FALSE]
      Ysub <- ydata[seq(1, nrow(xdata))[!seq(1, nrow(xdata)) %in% s], , drop = FALSE]
      mybeta2 <- SelectionAlgo(
        xdata = Xsub, ydata = Ysub,
        Lambda = Lambda[, 1], group_x = group_x,
        family = family, implementation = implementation, ...
      )
      
      # Resampling if model failed to converge
      while (is.infinite(mybeta1$selected[1]) | is.infinite(mybeta2$selected[1])) {
        s <- Resample(data = ydata, family = family, tau = tau, resampling = resampling, ...)
        
        # First subset
        Xsub <- xdata[s, , drop = FALSE]
        Ysub <- ydata[s, , drop = FALSE]
        mybeta <- SelectionAlgo(
          xdata = Xsub, ydata = Ysub,
          Lambda = Lambda[, 1], group_x = group_x,
          family = family, implementation = implementation, ...
        )
        
        # Complementary subset
        Xsub <- xdata[seq(1, nrow(xdata))[!seq(1, nrow(xdata)) %in% s], , drop = FALSE]
        Ysub <- ydata[seq(1, nrow(xdata))[!seq(1, nrow(xdata)) %in% s], , drop = FALSE]
        mybeta <- SelectionAlgo(
          xdata = Xsub, ydata = Ysub,
          Lambda = Lambda[, 1], group_x = group_x,
          family = family, implementation = implementation, ...
        )
      }
      
      # Storing beta coefficients from first set
      Beta[rownames(mybeta1$selected), colnames(mybeta1$selected), k] <- mybeta1$selected
      
      # Storing all beta coefficients from first set
      if (length(dim(Beta_full)) == 3) {
        Beta_full[rownames(mybeta1$beta_full), colnames(mybeta1$beta_full), k] <- mybeta1$beta_full
      } else {
        Beta_full[rownames(mybeta1$beta_full), colnames(mybeta1$beta_full), k, ] <- mybeta1$beta_full
      }
      
      # Storing beta coefficients from complementary set
      Beta[rownames(mybeta2$selected), colnames(mybeta2$selected), ceiling(K / 2) + k] <- mybeta2$selected
      
      # Storing all beta coefficients from complementary set
      if (length(dim(Beta_full)) == 3) {
        Beta_full[rownames(mybeta2$beta_full), colnames(mybeta2$beta_full), ceiling(K / 2) + k] <- mybeta2$beta_full
      } else {
        Beta_full[rownames(mybeta2$beta_full), colnames(mybeta2$beta_full), ceiling(K / 2) + k, ] <- mybeta2$beta_full
      }
      
      if (verbose) {
        utils::setTxtProgressBar(pb, 2 * k / K)
      }
    }
    
    # Computing the simultaneous selection proportions
    bigstab <- matrix(0, nrow = nrow(Beta), ncol = ncol(Beta))
    colnames(bigstab) <- colnames(Beta)
    rownames(bigstab) <- rownames(Beta)
    n_valid <- rep(0, nrow(Beta))
    for (k in 1:ceiling(K / 2)) {
      A1 <- ifelse(Beta[, , k] != 0, yes = 1, no = 0)
      A2 <- ifelse(Beta[, , ceiling(K / 2) + k] != 0, yes = 1, no = 0)
      A <- A1 + A2
      A <- ifelse(A == 2, yes = 1, no = 0)
      tmp_missing <- apply(A, 1, FUN = function(x) {
        any(is.na(x))
      })
      A[which(tmp_missing), ] <- 0
      n_valid <- n_valid + ifelse(!tmp_missing, yes = 1, no = 0)
      bigstab <- bigstab + A
    }
    bigstab <- bigstab / n_valid
  }
  
  if (verbose) {
    cat("\n")
  }
  
  # Computation of the stability score over Lambda and pi_list
  if (group_penalisation) {
    metrics <- StabilityMetrics(
      selprop = bigstab, pk = NULL, pi_list = pi_list, K = K, n_cat = n_cat,
      Sequential_template = NULL, graph = FALSE, group = group_x,
      PFER_method = PFER_method, PFER_thr_blocks = PFER_thr, FDP_thr_blocks = FDP_thr
    )
  } else {
    metrics <- StabilityMetrics(
      selprop = bigstab, pk = NULL, pi_list = pi_list, K = K, n_cat = n_cat,
      Sequential_template = NULL, graph = FALSE,
      PFER_method = PFER_method, PFER_thr_blocks = PFER_thr, FDP_thr_blocks = FDP_thr
    )
  }
  if (verbose) {
    utils::setTxtProgressBar(pb, 1)
    cat("\n")
  }
  Beta <- Beta_full
  
  # Preparing outputs
  myimplementation <- as.character(substitute(implementation, env = parent.frame(n = 2)))
  if (is.function(resampling)) {
    myresampling <- as.character(substitute(resampling))
  } else {
    myresampling <- resampling
  }
  out <- list(
    S = metrics$S, Lambda = Lambda,
    Q = metrics$Q, Q_s = metrics$Q_s, P = metrics$P,
    PFER = metrics$PFER, FDP = metrics$FDP,
    S_2d = metrics$S_2d, PFER_2d = metrics$PFER_2d, FDP_2d = metrics$FDP_2d,
    selprop = bigstab, Beta = Beta,
    methods = list(
      type = "variable_selection", implementation = myimplementation, family = family,
      resampling = myresampling, cpss = cpss, PFER_method = PFER_method
    ),
    params = list(
      K = K, pi_list = pi_list, tau = tau, n_cat = n_cat,
      pk = ncol(xdata), n = nrow(xdata),
      PFER_thr = PFER_thr, FDP_thr = FDP_thr,
      seed = seed
    )
  )
  
  if (output_data) {
    out$params <- c(out$params, list(xdata = xdata, ydata = ydata))
  }
  
  # Defining the class
  class(out) <- "variable_selection"
  
  return(out)
}

VariableSelection <- function(xdata, ydata = NULL, Lambda = NULL, pi_list = seq(0.01, 0.99, by = 0.01),
                              K = 100, tau = 0.5, seed = 1, n_cat = NULL,
                              family = "gaussian", implementation = PenalisedRegression,
                              resampling = "subsampling", cpss = FALSE,
                              PFER_method = "MB", PFER_thr = Inf, FDP_thr = Inf,
                              Lambda_cardinal = 100, group_x = NULL, group_penalisation = FALSE,
                              n_cores = 1, output_data = FALSE, verbose = TRUE, beep = NULL, ...) {
  # Defining Lambda if used with sparse PCA or PLS
  if (is.null(Lambda)) {
    if (as.character(substitute(implementation)) %in% c("SparseGroupPLS", "GroupPLS")) {
      Lambda <- seq(1, length(group_x) - 1)
    }
    if (as.character(substitute(implementation)) %in% c("SparsePLS", "SparsePCA")) {
      Lambda <- seq(1, ncol(xdata) - 1)
    }
  }
  
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
  
  # Checking that group_x is provided for group penalisation
  if (group_penalisation) {
    if (is.null(group_x)) {
      stop("Please provide argument 'group_x' for group penalisation. Argument 'group_x' should be a vector with the number of variables in each group.")
    }
  }
  
  if (is.null(Lambda)) {
    # Defining grid of lambda values (using glmnet implementation)
    Lambda <- LambdaGridRegression(
      xdata = xdata, ydata = ydata, tau = tau, seed = seed,
      family = family,
      resampling = resampling,
      Lambda_cardinal = Lambda_cardinal, check_input = FALSE, ...
    )
  }
  
  # Check if parallelisation is possible (forking)
  if (.Platform$OS.type != "unix") {
    if (n_cores > 1) {
      warning("Invalid input for argument 'n_cores'. Parallelisation relies on forking, it is only available on Unix systems.")
    }
    n_cores <- 1
  }
  
  # Stability selection and score
  mypar <- parallel::mclapply(X = 1:n_cores, FUN = function(k) {
    return(SerialRegression(
      xdata = xdata, ydata = ydata, Lambda = Lambda, pi_list = pi_list,
      K = ceiling(K / n_cores), tau = tau, seed = as.numeric(paste0(seed, k)), n_cat = n_cat,
      family = family, implementation = implementation,
      resampling = resampling, cpss = cpss,
      PFER_method = PFER_method, PFER_thr = PFER_thr, FDP_thr = FDP_thr,
      group_x = group_x, group_penalisation = group_penalisation,
      output_data = output_data, verbose = verbose, ...
    ))
  })
  
  # Combining the outputs from parallel iterations
  out <- mypar[[1]]
  if (n_cores > 1) {
    for (i in 2:length(mypar)) {
      out <- do.call(Combine, list(stability1 = out, stability2 = mypar[[i]]))
    }
  }
  
  # Re-set the function names
  if ("methods" %in% names(out)) {
    myimplementation <- as.character(substitute(implementation))
    if (is.function(resampling)) {
      myresampling <- as.character(substitute(resampling))
    } else {
      myresampling <- resampling
    }
    out$methods$implementation <- myimplementation
    out$methods$resampling <- myresampling
  }
  
  # Defining the class
  class(out) <- "variable_selection"
  
  # Making beep
  if (!is.null(beep)) {
    beepr::beep(sound = beep)
  }
  
  return(out)
}

CheckParamRegression <- function(Lambda = NULL, pi_list = seq(0.6, 0.9, by = 0.01),
                                 K = 100, tau = 0.5, seed = 1, n_cat = NULL,
                                 family = "gaussian", implementation = PenalisedRegression,
                                 resampling = "subsampling", PFER_method = "MB", PFER_thr = Inf, FDP_thr = Inf,
                                 Lambda_cardinal = 100,
                                 verbose = TRUE) {
  # List of arguments
  myargs <- c(
    "Lambda", "pi_list", "K", "tau", "seed", "n_cat",
    "family",
    "PFER_method", "PFER_thr", "FDP_thr",
    "Lambda_cardinal", "verbose"
  )
  
  # Checking the inputs (Lambda)
  if (!is.null(Lambda)) {
    if (is.matrix(Lambda)) {
      Lambda_copy <- Lambda
      Lambda <- NULL
      for (k in 1:ncol(Lambda_copy)) {
        Lambda <- cbind(Lambda, as.numeric(Lambda_copy[, k]))
      }
    } else {
      Lambda <- as.numeric(Lambda)
      Lambda <- cbind(Lambda)
    }
    if (any(is.na(Lambda))) {
      if (all(is.na(Lambda))) {
        stop("Invalid input for argument 'Lambda'. The input only contains missing values.")
      } else {
        Lambda <- as.matrix(stats::na.exclude(Lambda))
        warning("Invalid input for argument 'Lambda'. The input contains missing values. These have been excluded.")
      }
    }
    rownames(Lambda) <- paste0("s", seq(0, nrow(Lambda) - 1))
  }
  
  # Checking the inputs (pi_list)
  pi_list <- sort(pi_list)
  restricted <- FALSE
  if (!is.null(n_cat)) {
    if (n_cat == 3) {
      restricted <- TRUE
    }
  }
  if (restricted) {
    if (any(pi_list > 0.5) & any(pi_list < 1)) {
      if ((min(pi_list) < 0.5) | (max(pi_list) > 1)) {
        warning("The values in 'pi_list' must be between 0.5 and 1. All other values were discarded.")
        pi_list <- pi_list[which((pi_list > 0.5) & (pi_list < 1))]
      }
    } else {
      stop("Invalid input for argument 'pi_list'. The values in the vector must be greater than 0.5 and lower than 1. To consider thresholds below 0.5, argument 'n_cat' must be set to 2.")
    }
  } else {
    if (any(pi_list > 0) & any(pi_list < 1)) {
      if ((min(pi_list) < 0) | (max(pi_list) > 1)) {
        warning("The values in 'pi_list' must be between 0 and 1. All other values were discarded.")
        pi_list <- pi_list[which((pi_list > 0) & (pi_list < 1))]
      }
    } else {
      stop("Invalid input for argument 'pi_list'. The values in the vector must be greater than 0 and lower than 1.")
    }
  }
  
  # Checking the inputs (K)
  K <- as.numeric(K)
  if ((length(K) != 1) | is.na(K)) {
    warning("Invalid input for argument 'K'. The number of resampling iterations 'K' must be a single number.")
    K <- 100
  }
  
  # Checking the inputs (tau)
  tau <- as.numeric(tau)
  if ((length(tau) != 1) | is.na(tau) | (tau >= 1) | (tau <= 0)) {
    warning("Invalid input for argument 'tau'. The subsample size 'tau' must be a number between 0 and 1. The default value (0.5) was used.")
    tau <- 0.5
  }
  
  # Checking the inputs (seed)
  seed <- as.numeric(seed)
  if ((length(seed) != 1) | is.na(seed)) {
    warning("Invalid input for argument 'seed'. The argument 'seed' must be a single number. The default value (1) was used.")
    seed <- 1
  }
  
  # Checking the inputs (n_cat)
  if (!is.null(n_cat)) {
    n_cat <- as.numeric(n_cat)
    if ((length(n_cat) != 1) | is.na(n_cat)) {
      warning("Invalid input for argument 'n_cat'. The argument 'n_cat' must be set to 2 or 3. The default value (3) was used.")
      n_cat <- 3
    }
  }
  
  # Checking the inputs (family)
  family <- as.character(family)
  if ((length(family) != 1) | is.na(family)) {
    stop("Invalid input for argument 'family'. The argument 'family' must be a character string.")
  }
  
  # Checking the inputs (implementation)
  if (!is.function(implementation)) {
    stop("Invalid input for argument 'implementation'. This argument must be a function to use for variable selection.")
  }
  
  # Checking the inputs (resampling)
  if ((!is.function(resampling)) & (!is.character(resampling))) {
    stop("Invalid input for argument 'resampling'. The argument 'resampling' must be a character string. Possible values are: 'subsampling', 'bootstrap' or the name of a function.")
  }
  
  # Checking the inputs (PFER_method)
  PFER_method <- as.character(PFER_method)
  if ((length(PFER_method) != 1) | (!PFER_method %in% c("MB", "SS"))) {
    stop("Invalid input for argument 'PFER_method'. Possible values are: 'MB' or 'SS'.")
  }
  
  # Checking the inputs (PFER_method and resampling)
  if (is.character(resampling)) {
    if ((PFER_method == "SS") & (resampling == "bootstrap")) {
      warning("Arguments 'resampling' and 'PFER_method' are not compatible. With 'PFER_method' set to 'SS', the resampling is done with complementary pairs of subsamples.")
      resampling <- "subsampling"
    }
  }
  
  # Checking the inputs (PFER_thr)
  PFER_thr <- as.numeric(PFER_thr)
  if ((length(PFER_thr) != 1) | any(is.na(PFER_thr)) | any(PFER_thr <= 0)) {
    warning("Invalid input for argument 'PFER_thr'. The threshold in the upper-bound of the expected number of False Positives 'PFER_thr' must be a single positive number (or Inf). The default value (Inf) was used.")
    PFER_thr <- Inf
  }
  
  # Checking the inputs (FDP_thr)
  FDP_thr <- as.numeric(FDP_thr)
  if ((length(FDP_thr) != 1) | any(is.na(FDP_thr)) | any((!is.infinite(FDP_thr)) & (FDP_thr <= 0)) | any((!is.infinite(FDP_thr)) & (FDP_thr > 1))) {
    warning("Invalid input for argument 'FDP_thr'. The threshold in the upper-bound of the False Discovery Proportion 'FDP_thr' must be a single number between 0 and 1 (or Inf to deactivate). The default value (Inf) was used.")
    FDP_thr <- Inf
  }
  
  # Checking the inputs (PFER_thr and FDP_thr)
  if ((!is.infinite(PFER_thr)) & (!is.infinite(FDP_thr))) {
    warning("Arguments 'PFER_thr' and 'FDP_thr' are not compatible. Only one of these two arguments can be used (i.e. not set to Inf). Argument 'PFER_thr' was used.")
    FDP_thr <- Inf
  }
  
  # Checking the inputs (Lambda_cardinal)
  Lambda_cardinal <- as.numeric(Lambda_cardinal)
  if (is.null(Lambda)) {
    if ((length(Lambda_cardinal) != 1) | is.na(Lambda_cardinal) | (Lambda_cardinal < 1)) {
      warning("Invalid input for argument 'Lambda_cardinal'. The argument 'Lambda_cardinal' must be a single positive number. A value of 10 was used.")
      Lambda_cardinal <- 10
    }
  }
  
  # Checking the inputs (verbose)
  verbose <- as.logical(verbose)
  if ((length(verbose) != 1) | is.na(verbose)) {
    warning("Invalid input for argument 'verbose'. The argument 'verbose' must be logical (TRUE or FALSE). The default value (TRUE) was used.")
    verbose <- TRUE
  }
  
  # Assigning checked values to the parent function
  for (i in 1:length(myargs)) {
    if (!is.null(get(myargs[i]))) {
      assign(myargs[i], get(myargs[i]), envir = parent.frame(n = 1))
    }
  }
  
}


CheckDataRegression <- function(xdata, ydata = NULL,
                                family = "gaussian",
                                verbose = TRUE) {
  # List of arguments
  myargs <- c(
    "xdata", "ydata", "family"
  )
  
  # Checking the inputs (xdata and ydata)
  xdata <- as.matrix(xdata)
  if (!is.null(ydata)) {
    if (sum(is.na(xdata)) > 0) {
      stop("Invalid input for argument 'xdata'. Missing values are not allowed in 'xdata'.")
    }
    if (sum(is.na(ydata)) > 0) {
      stop("Invalid input for argument 'ydata'. Missing values are not allowed in 'ydata'.")
    }
    if (nrow(xdata) < 10) {
      stop("Invalid input for argument 'xdata'. Not enough data.")
    }
  }
  
  # Preparing xdata
  if (is.null(colnames(xdata))) {
    colnames(xdata) <- paste0("var", 1:ncol(xdata))
  }
  
  # Preparing ydata
  if (!is.null(ydata)) {
    # Turning data frame into matrix
    if (is.data.frame(ydata)) {
      ydata <- as.matrix(ydata)
    }
    # Turning matrix with one column into vector
    if (is.matrix(ydata)) {
      if (ncol(ydata) == 1) {
        ydata <- as.vector(ydata)
      }
    }
    # Turning vector into factor
    if ((family %in% c("binomial", "multinomial")) & is.vector(ydata)) {
      ydata <- as.factor(ydata)
    }
    
    # Defining reference category and final data type
    if (is.factor(ydata)) {
      if ((family %in% c("binomial", "multinomial")) & verbose) {
        message(paste0("Reference category: ", levels(ydata)[1]))
        message(paste0("Other categorie(s): ", paste(levels(ydata)[-1], collapse = ", ")))
      }
      ydata <- as.numeric(ydata) - 1
    }
    ydata <- matrix(ydata, ncol = 1)
  }
  
  # Checking the inputs (xdata and ydata)
  if (!is.null(ydata)) {
    if (nrow(xdata) != nrow(ydata)) {
      stop("Arguments 'xdata' and 'ydata' are not compatible. They have different numbers of observations.")
    }
  }
  
  # Creating dummy ydata (for resampling in unsupervised models)
  if (is.null(ydata)) {
    ydata <- cbind(rep(0, nrow(xdata)))
  }
  
  # Naming rows of xdata and ydata
  if (is.null(rownames(xdata)) & is.null(rownames(ydata))) {
    rownames(xdata) <- paste0("obs", 1:nrow(xdata))
    rownames(ydata) <- rownames(xdata)
  } else {
    if ((is.null(rownames(xdata))) & (!is.null(rownames(ydata)))) {
      rownames(xdata) <- rownames(ydata)
    }
    if ((!is.null(rownames(xdata))) & (is.null(rownames(ydata)))) {
      rownames(ydata) <- rownames(xdata)
    }
  }
  
  # Re-ordering the datasets to ensure that subsamples will be the same regardless of the order of observations in the input
  ids <- sort.list(rownames(xdata))
  xdata <- xdata[ids, , drop = FALSE]
  ydata <- ydata[ids, , drop = FALSE]
  
  # Further checking/preparing ydata
  if ((family == "cox")) {
    if ((ncol(ydata) != 2) | (length(unique(ydata[, 2])) != 2)) {
      stop("Invalid input for argument 'ydata'. For Cox regression using glmnet, the argument 'ydata' needs to be a matrix or data frame with two columns: the time to event and binary status.")
    }
    colnames(ydata) <- c("time", "status")
    tmp <- as.factor(ydata[, 2])
    if (verbose) {
      message(paste0("Reference category: ", levels(tmp)[1]))
      message(paste0("Other category: ", levels(tmp)[2]))
    }
    ydata[, 2] <- as.numeric(tmp) - 1
    ydata <- as.matrix(ydata)
  }
  if ((family %in% c("binomial", "multinomial"))) {
    if (ncol(ydata) > 1) {
      ydata <- DummyToCategories(x = ydata, verbose = verbose)
      ydata <- matrix(ydata, ncol = 1)
      rownames(ydata) <- rownames(xdata)
    }
    ytmp <- as.numeric(table(ydata))
    if (any(ytmp == 1)) {
      stop("At least one category in 'ydata' with only one observation.")
    }
  }
  
  # Assigning checked values to the parent function
  for (i in 1:length(myargs)) {
    if (!is.null(get(myargs[i]))) {
      assign(myargs[i], get(myargs[i]), envir = parent.frame(n = 1))
    }
  }
}
 
  Heatmap <- function(mat, col = c("ivory", "navajowhite", "tomato", "darkred"),
                     resolution = 10000, bty = "o",
                     axes = TRUE, cex.axis = 1, xlas = 2, ylas = 2,
                     text = FALSE, cex = 1,
                     legend = TRUE, legend_length = NULL, legend_range = NULL, cex.legend = 1, ...) {
   oldpar <- graphics::par("xpd", "xaxs", "yaxs", no.readonly = TRUE)
   on.exit(graphics::par(oldpar))
   
   # Transposing the input matrix so that rows are rows
   mat <- t(mat)
   
   # Defining the legend length
   if (is.null(legend_length)) {
     legend_length <- ncol(mat)
   }
   
   # Preparing colours
   col <- grDevices::colorRampPalette(col)(resolution)
   names(col) <- 1:resolution
   
   # Re-formatting matrix
   mat <- mat[, ncol(mat):1, drop = FALSE]
   vect <- as.vector(mat)
   
   # Defining extreme values
   if (is.null(legend_range)) {
     # myrange <- c(min(vect, na.rm = TRUE), max(vect, na.rm = TRUE))
     myrange <- range(vect, na.rm = TRUE)
     myrange <- c(floor(myrange[1]), ceiling(myrange[2]))
   } else {
     myrange <- legend_range
   }
   
   # Getting corresponding colours
   mycol <- as.character(cut(vect, breaks = seq(myrange[1], myrange[2], length.out = resolution + 1), labels = 1:resolution, include.lowest = TRUE))
   mycol_mat <- matrix(mycol, ncol = ncol(mat))
   
   # Making heatmap
   withr::local_par(xaxs = "i", yaxs = "i")
   plot(NA,
        xlim = c(0, nrow(mycol_mat)), ylim = c(0, ncol(mycol_mat)),
        xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n"
   )
   for (i in 0:(nrow(mycol_mat) - 1)) {
     for (j in 0:(ncol(mycol_mat) - 1)) {
       graphics::polygon(
         x = c(i, i + 1, i + 1, i), y = c(j, j, j + 1, j + 1),
         col = col[mycol_mat[i + 1, j + 1]],
         border = col[mycol_mat[i + 1, j + 1]]
       )
     }
   }
   if (axes) {
     if (!is.null(rownames(mat))) {
       graphics::axis(side = 1, at = 1:nrow(mat) - 0.5, labels = rownames(mat), las = xlas, cex.axis = cex.axis)
     }
     if (!is.null(colnames(mat))) {
       graphics::axis(side = 2, at = 1:ncol(mat) - 0.5, labels = colnames(mat), las = ylas, cex.axis = cex.axis)
     }
   }
   if (bty == "o") {
     graphics::box()
   }
   
   # Showing text
   if (text) {
     for (i in 1:nrow(mat)) {
       for (j in 1:ncol(mat)) {
         text(i - 0.5, j - 0.5,
              cex = cex,
              labels = formatC(mat[i, j], ...)
         )
       }
     }
   }
   
   # Adding colour bar (legend)
   if (legend) {
     withr::local_par(list(xpd = TRUE))
     legend_width_factor <- 1.05
     mylegend_values <- grDevices::axisTicks(c(myrange[1], myrange[2]), log = FALSE)
     mylegend_ids <- as.numeric(as.character(cut(mylegend_values,
                                                 breaks = seq(myrange[1], myrange[2], length.out = resolution + 1),
                                                 labels = 1:resolution, include.lowest = TRUE
     )))
     ypos <- ncol(mat)
     xpos <- nrow(mat) * 1.05
     for (l in 1:length(col)) {
       graphics::polygon(
         x = c(xpos, xpos * legend_width_factor, xpos * legend_width_factor, xpos),
         y = c(
           ypos - legend_length + legend_length * l / length(col),
           ypos - legend_length + legend_length * l / length(col),
           ypos - legend_length + legend_length * (l + 1) / length(col),
           ypos - legend_length + legend_length * (l + 1) / length(col)
         ),
         col = col[l], border = col[l]
       )
       if (l %in% mylegend_ids) {
         graphics::text(
           x = xpos * legend_width_factor, y = ypos - legend_length + legend_length * (l + 0.5) / length(col),
           labels = paste0("- ", mylegend_values[which(mylegend_ids == l)]), adj = c(0, 0.5), cex = cex.legend
         )
       }
     }
     withr::local_par(list(xpd = FALSE)) # for legend
   }
 }
 
 CalibrationPlot <- function(stability, block_id = NULL,
                             col = NULL,
                             pch = 19, cex = 0.7,
                             xlim = NULL, ylim = NULL, bty = "o",
                             lines = TRUE, lty = 3, lwd = 2,
                             show_argmax = TRUE,
                             show_pix = FALSE, show_piy = FALSE, offset = 0.3,
                             legend = TRUE, legend_length = NULL, legend_range = NULL, ncol = 1,
                             xlab = NULL, ylab = NULL, zlab = expression(italic(q)),
                             xlas = 2, ylas = NULL, zlas = 2,
                             cex.lab = 1.5, cex.axis = 1, cex.legend = 1.2,
                             xgrid = FALSE, ygrid = FALSE,
                             params = c("ny", "alphay", "nx", "alphax")) {
   oldpar <- graphics::par("xpd", "xaxs", "yaxs", no.readonly = TRUE)
   on.exit(graphics::par(oldpar))
   
   # To deal with later: showing calibration of clustering or selection
   # clustering <- FALSE
   clustering <- ifelse(inherits(stability, "clustering"), yes = TRUE, no = FALSE)
   heatmap <- TRUE
   
   if (clustering) {
     ylas <- 1
     CalibrationCurve(
       stability = stability, bty = bty, xlab = xlab, ylab = ylab,
       cex.lab = cex.lab, cex.axis = cex.axis, cex.legend = cex.legend,
       pch = pch, lines = lines, col = col, legend = legend, ncol = ncol
     )
   } else {
     ylas <- 0
     
     if (inherits(stability, "bi_selection")) {
       # Extracting summary information
       x <- stability$summary_full
       
       # Checking input
       params <- unique(params)
       all_params <- colnames(stability$summary)
       all_params <- all_params[!all_params %in% c("comp", "S", "pix", "piy")]
       if (any(!all_params %in% params)) {
         params <- all_params
         warning(paste0(
           "Invalid input for argument 'params'. Please provide a vector with all the following: ",
           paste(all_params, collapse = ", "), "."
         ))
       }
       params <- params[params %in% all_params]
       
       # Identifying parameters
       params <- params[params %in% colnames(x)]
       
       # Defining default arguments
       if (is.null(ylab)) {
         ylab <- "Stability Score"
       }
       
       if (is.null(xlab)) {
         if (length(params) > 1) {
           xlab <- ""
         } else {
           xlab <- expression(n[X])
         }
       }
       
       if (is.null(col)) {
         col <- grDevices::colorRampPalette(c("navy", "darkred"))(nrow(stability$summary))
       } else {
         col <- grDevices::colorRampPalette(col)(nrow(stability$summary))
       }
       
       if (length(unique(x$comp)) == 1) {
         legend <- FALSE
       }
       
       if (is.null(xlim)) {
         xlim <- c(0.5, max(sapply(split(x, f = x$comp), nrow)) + 0.5)
       }
       
       if (is.null(ylim)) {
         ylim <- range(x$S, na.rm = TRUE)
         if (legend) {
           ylim[2] <- ylim[2] + diff(ylim) * 0.15
         }
       }
       
       # Drawing one set of points per component
       for (comp_id in unique(x$comp)) {
         tmp <- x[which(x$comp == comp_id), ]
         
         # Ensuring increasing ny
         tmp <- tmp[do.call(order, tmp[, params, drop = FALSE]), ]
         # if ("ny" %in% colnames(tmp)) {
         #   tmp=tmp[order(tmp$ny, tmp$nx), ]
         # }
         # tmp=tmp[order(lapply(params, FUN=function(param_id){with(tmp, eval(parse(text=param_id)))})),]
         # if ("alphax" %in% colnames(tmp))
         
         if (comp_id == min(x$comp)) {
           # Initialising the plot
           plot(NA,
                xlim = xlim, ylim = ylim, bty = bty,
                xlab = xlab, ylab = ylab, cex.lab = cex.lab,
                cex.axis = cex.axis,
                xaxt = "n", las = ylas
           )
           
           # Defining vertical grid
           if (xgrid) {
             withr::local_par(list(xpd = FALSE))
             graphics::abline(v = 1:nrow(tmp), lty = 3, col = "grey")
           }
           
           # Defining horizontal grid
           if (ygrid) {
             withr::local_par(list(xpd = FALSE))
             graphics::abline(h = graphics::axTicks(side = 2), lty = 3, col = "grey")
           }
           
           # Adding x-axes
           for (param_id in 1:length(params)) {
             if (param_id == 1) {
               graphics::axis(
                 side = 1, at = 1:nrow(tmp),
                 labels = tmp[, rev(params)[param_id]],
                 cex.axis = cex.axis, las = xlas
               )
             } else {
               ids <- c(1, which(diff(tmp[, rev(params)[param_id]]) != 0) + 1)
               ids <- c(ids - 0.5, nrow(tmp) + 0.5)
               graphics::axis(side = 1, at = ids, labels = NA, line = (param_id - 1) * 3)
               withr::local_par(list(xpd = FALSE))
               graphics::abline(v = ids, lty = 2)
               ids <- apply(rbind(ids[-1], ids[-length(ids)]), 2, mean)
               graphics::axis(
                 side = 1, at = ids, labels = tmp[ids, rev(params)[param_id]],
                 line = (param_id - 1) * 3, tick = FALSE, cex.axis = cex.axis, las = xlas
               )
             }
           }
           # graphics::axis(side = 1, at = 1:nrow(tmp), labels = tmp$nx, cex.axis = cex.axis, las = xlas)
           # if ("ny" %in% colnames(tmp)) {
           #   ids <- c(which(!duplicated(tmp$ny)) - 0.5, nrow(tmp) + 0.5)
           #   graphics::axis(side = 1, at = ids, labels = NA, line = 3)
           #   withr::local_par(list(xpd = FALSE))
           #   graphics::abline(v = ids, lty = 2)
           #   ids <- apply(rbind(ids[-1], ids[-length(ids)]), 2, mean)
           #   graphics::axis(side = 1, at = ids, labels = unique(tmp$ny), line = 3, tick = FALSE, cex.axis = cex.axis, las = xlas)
           # }
           
           # Adding x-labels
           if (length(params) > 1) {
             for (param_id in 1:length(params)) {
               if (rev(params)[param_id] == "nx") {
                 graphics::mtext(text = expression(n[X]), side = 1, at = 0.5 - nrow(tmp) * 0.1, line = (param_id - 1) * 3 + 1, cex = cex.lab)
               }
               if (rev(params)[param_id] == "alphax") {
                 graphics::mtext(text = expression(alpha[X]), side = 1, at = 0.5 - nrow(tmp) * 0.1, line = (param_id - 1) * 3 + 1, cex = cex.lab)
               }
               if (rev(params)[param_id] == "ny") {
                 graphics::mtext(text = expression(n[Y]), side = 1, at = 0.5 - nrow(tmp) * 0.1, line = (param_id - 1) * 3 + 1, cex = cex.lab)
               }
               if (rev(params)[param_id] == "alphay") {
                 graphics::mtext(text = expression(alpha[Y]), side = 1, at = 0.5 - nrow(tmp) * 0.1, line = (param_id - 1) * 3 + 1, cex = cex.lab)
               }
             }
             # graphics::mtext(text = expression(n[X]), side = 1, at = -nrow(tmp) * 0.06, line = 1, cex = cex.lab)
             # graphics::mtext(text = expression(n[Y]), side = 1, at = -nrow(tmp) * 0.06, line = 4, cex = cex.lab)
           }
         }
         
         # Adding calibrated lines
         if (show_argmax) {
           withr::local_par(list(xpd = FALSE))
           graphics::abline(v = which.max(tmp$S), lty = 3, col = col[comp_id])
         }
         
         # Adding lines
         if (lines) {
           # if ("ny" %in% colnames(tmp)) {
           #   for (y_value in unique(tmp$ny)) {
           #     graphics::lines(which(tmp$ny == y_value),
           #                     tmp[which(tmp$ny == y_value), "S"],
           #                     col = col[comp_id],
           #                     lty = lty, lwd = lwd
           #     )
           #   }
           # } else {
           graphics::lines(1:nrow(tmp),
                           tmp$S,
                           col = col[comp_id],
                           lty = lty, lwd = lwd
           )
           # }
         }
         
         # Adding data points
         graphics::points(tmp$S,
                          pch = pch,
                          col = col[comp_id],
                          cex = cex
         )
         
         # Adding pi values
         if ((show_pix) & (!show_piy)) {
           graphics::text(tmp$S,
                          labels = tmp$pix,
                          col = col[comp_id],
                          cex = cex, pos = 3,
                          offset = offset
           )
         }
         
         if ((!show_pix) & (show_piy)) {
           graphics::text(tmp$S,
                          labels = tmp$piy,
                          col = col[comp_id],
                          cex = cex, pos = 3,
                          offset = offset
           )
         }
         
         if ((show_pix) & (show_piy)) {
           for (k in 1:nrow(tmp)) {
             graphics::text(k, tmp[k, "S"],
                            labels = eval(parse(text = paste0("expression(pi[x]*' = ", tmp[k, "pix"], " ; '*pi[y]*' = ", tmp[k, "piy"], "')"))),
                            col = col[comp_id],
                            cex = cex, pos = 3,
                            offset = offset
             )
           }
         }
       }
       
       # Adding legend
       if (legend) {
         graphics::legend("top",
                          col = col, lty = lty, pch = pch, lwd = lwd,
                          legend = paste0("Component ", unique(x$comp)),
                          horiz = TRUE, bg = "white", cex = cex.legend
         )
       }
     } else {
       # Defining default arguments
       if (heatmap) {
         metric <- "both"
         if (is.null(col)) {
           col <- c("ivory", "navajowhite", "tomato", "darkred")
         }
         if (is.null(ylab)) {
           if (clustering) {
             ylab <- expression(n[c])
           } else {
             ylab <- expression(pi)
           }
         }
       } else {
         metric <- "lambda"
         if (is.null(col)) {
           col <- "navy"
         }
         if (is.null(ylab)) {
           ylab <- "Stability Score"
         }
       }
       if (is.null(xlab)) {
         xlab <- expression(lambda)
       }
       
       # Extracting the number of blocks/components
       if ((stability$methods$type == "graphical_model") & (is.null(block_id))) {
         bigblocks <- BlockMatrix(stability$params$pk)
         nblocks <- length(stability$params$pk) * (length(stability$params$pk) + 1) / 2
         bigblocks_vect <- factor(bigblocks[upper.tri(bigblocks)], levels = 1:nblocks)
         N_blocks <- unname(table(bigblocks_vect))
         blocks <- levels(bigblocks_vect)
         names(N_blocks) <- blocks
         block_id <- which(N_blocks != 0)
       } else {
         if (is.null(block_id)) {
           block_id <- 1
         }
       }
       nblocks <- length(block_id)
       
       if (metric == "both") {
         for (b in block_id) {
           # Printing block ID
           if (length(block_id) > 1) {
             message(paste0("Block ", b))
           }
           
           # Extracting the stability scores
           if (clustering) {
             mat <- matrix(stability$Sc, ncol = length(unique(stability$Lambda)))
             rownames(mat) <- formatC(unique(stability$nc), format = "f", digits = 0)
             colnames(mat) <- formatC(unique(stability$Lambda), format = "e", digits = 2)
             mat <- t(mat)
             mat <- mat[nrow(mat):1, ncol(mat):1]
             ids <- which(apply(mat, 1, FUN = function(x) {
               any(!is.na(x))
             }))
             mat <- mat[ids, , drop = FALSE]
           } else {
             if (length(stability$params$pk) == 1) {
               mat <- stability$S_2d
               ids <- which(apply(mat, 1, FUN = function(x) {
                 any(!is.na(x))
               }))
               mat <- mat[ids, , drop = FALSE]
             } else {
               mat <- stability$S_2d[, , b]
               ids <- which(apply(mat, 1, FUN = function(x) {
                 any(!is.na(x))
               }))
               mat <- mat[ids, , drop = FALSE]
             }
             
             # Setting row and column names
             colnames(mat) <- stability$params$pi_list
             if (grepl("penalised", tolower(stability$methods$implementation))) {
               rownames(mat) <- formatC(stability$Lambda[, b], format = "e", digits = 2)[ids]
             } else {
               rownames(mat) <- (stability$Lambda[, b])[ids]
             }
           }
           
           # Extracting corresponding numbers of selected variables (q)
           Q <- stability$Q[, b]
           Q <- Q[ids]
           
           # Heatmap representation
           Heatmap(t(mat[nrow(mat):1, ncol(mat):1]),
                   col = col, bty = bty, axes = FALSE,
                   text = FALSE, cex = 1, digits = 2,
                   legend = legend, legend_length = legend_length, legend_range = legend_range
           )
           
           # Adding calibrated lines
           if (show_argmax) {
             withr::local_par(list(xpd = FALSE))
             if (clustering) {
               myargmax <- Argmax(stability)
               print(myargmax)
               argmax_id <- c(
                 which(rownames(t(mat[nrow(mat):1, ncol(mat):1])) == formatC(myargmax[1], format = "f", digits = 0)),
                 which(colnames(t(mat[nrow(mat):1, ncol(mat):1])) == formatC(myargmax[2], format = "e", digits = 2))
               )
               print(argmax_id)
               graphics::abline(h = ncol(mat) - argmax_id[1] + 0.5, lty = 3)
               graphics::abline(v = argmax_id[2] - 0.5, lty = 3)
             } else {
               graphics::abline(v = nrow(mat) - which(stability$Lambda[ids, b] == Argmax(stability)[b, 1]) + 0.5, lty = 3)
               graphics::abline(h = which.min(abs(as.numeric(colnames(mat)) - Argmax(stability)[b, 2])) - 0.5, lty = 3)
             }
           }
           
           # Including axes
           if (clustering) {
             graphics::axis(
               side = 2, at = (1:ncol(mat)) - 0.5, las = ylas, cex.axis = cex.axis,
               labels = formatC(as.numeric(colnames(mat)), format = "f", digits = 0)
             )
           } else {
             graphics::axis(
               side = 2, at = (1:ncol(mat)) - 0.5, las = ylas, cex.axis = cex.axis,
               labels = formatC(as.numeric(colnames(mat)), format = "f", digits = 2)
             )
           }
           if (grepl("penalised", tolower(stability$methods$implementation))) {
             graphics::axis(
               side = 3, at = (1:nrow(mat)) - 0.5, las = zlas, cex.axis = cex.axis,
               labels = rev(formatC(Q, format = "f", big.mark = ",", digits = 0))
             )
             graphics::axis(side = 1, at = (1:nrow(mat)) - 0.5, las = xlas, labels = rev(rownames(mat)), cex.axis = cex.axis)
           } else {
             graphics::axis(side = 1, at = (1:nrow(mat)) - 0.5, las = xlas, labels = rev(rownames(mat)), cex.axis = cex.axis)
           }
           
           # Including axis labels
           graphics::mtext(text = ylab, side = 2, line = 3.5, cex = cex.lab)
           graphics::mtext(text = xlab, side = 1, line = 5.2, cex = cex.lab)
           if (grepl("penalised", tolower(stability$methods$implementation))) {
             graphics::mtext(text = zlab, side = 3, line = 3.5, cex = cex.lab)
           }
         }
       } else {
         if (metric == "lambda") {
           for (b in block_id) {
             # Extracting the stability scores
             if (length(stability$params$pk) == 1) {
               mat <- stability$S_2d
               ids <- which(apply(mat, 1, FUN = function(x) {
                 any(!is.na(x))
               }))
               mat <- mat[ids, , drop = FALSE]
             } else {
               mat <- stability$S_2d[, , b]
               ids <- which(apply(mat, 1, FUN = function(x) {
                 any(!is.na(x))
               }))
               mat <- mat[ids, , drop = FALSE]
             }
             
             # Extracting the best stability score (with optimal pi) for each lambda value
             vect <- apply(mat, 1, max, na.rm = TRUE)
             
             # Extracting corresponding numbers of selected variables (q)
             Q <- stability$Q[, b, drop = FALSE]
             Q <- Q[ids]
             
             # Extracting corresponding lambda values
             Lambda <- stability$Lambda[ids, b, drop = FALSE]
             
             # Re-ordering by decreasing lambda
             ids <- sort.list(Lambda, decreasing = TRUE)
             Lambda <- Lambda[ids]
             Q <- Q[ids]
             vect <- vect[ids]
             
             if (is.null(xlim)) {
               xlim <- range(Lambda, na.rm = TRUE)
             }
             
             if (is.null(ylim)) {
               ylim <- range(vect)
             }
             
             # Initialising the plot
             plot(NA,
                  xlim = xlim, ylim = ylim, bty = bty,
                  xlab = "", ylab = ylab, cex.lab = cex.lab,
                  cex.axis = cex.axis,
                  xaxt = "n", las = ylas
             )
             
             # Defining horizontal grid
             if (ygrid) {
               withr::local_par(list(xpd = FALSE))
               graphics::abline(h = graphics::axTicks(side = 2), lty = 3, col = "grey")
             }
             
             # Adding calibrated lines
             if (show_argmax) {
               withr::local_par(list(xpd = FALSE))
               graphics::abline(h = max(vect), lty = 3, col = col[1])
               graphics::abline(v = Lambda[which.max(vect)], lty = 3, col = col[1])
             }
             
             # Adding lines
             if (lines) {
               graphics::lines(Lambda, vect, col = col[1], lty = lty, lwd = lwd)
             }
             
             # Adding data points
             graphics::points(Lambda, vect, pch = pch, col = col[1], cex = cex)
             
             # Adding x-axis and z-axis and their labels
             lseq <- grDevices::axisTicks(range(Lambda, na.rm = TRUE), log = FALSE)
             xseq <- 1
             for (i in 1:length(lseq)) {
               xseq <- c(xseq, which.min(abs(Lambda - lseq[i])))
             }
             xseq <- c(xseq, length(Lambda))
             xseq <- unique(xseq)
             if (xgrid) {
               withr::local_par(list(xpd = FALSE))
               graphics::abline(v = Lambda[xseq], lty = 3, col = "grey")
             }
             graphics::axis(side = 1, at = Lambda[xseq], labels = formatC(Lambda[xseq], format = "e", digits = 2), las = xlas, cex.axis = cex.axis)
             graphics::axis(
               side = 3, at = Lambda[xseq], las = xlas,
               labels = formatC(Q[xseq], format = "f", big.mark = ",", digits = 0), cex.axis = cex.axis
             )
             graphics::mtext(text = xlab, side = 1, line = 5.2, cex = cex.lab)
             graphics::mtext(text = zlab, side = 3, line = 3.5, cex = cex.lab)
           }
         }
       }
     }
   }
 }
 


