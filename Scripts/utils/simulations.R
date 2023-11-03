#' Intermediate repository for utility functions to be implemented in fake
#' 
#' Data simulation for multivariate linear regression with pairwise interactions
#' 
#' Simulates data with outcome(s) and predictors, where only a subset of the
#' predictors and their multiplicative pairwise interactions
#' actually contributes to the definition of the outcome(s).
#' 
#' @param n number of observations in the simulated dataset. Not used if
#'   \code{xdata} is provided.
#' @param pk number of predictor variables. A subset of these variables
#'   contribute to the outcome definition (see argument \code{nu_xy}). Not used
#'   if \code{xdata} is provided.
#' @param xdata optional data matrix for the predictors with variables as
#'   columns and observations as rows. A subset of these variables contribute to
#'   the outcome definition (see argument \code{nu_xy}).
#' @param family type of regression model. Possible values include
#'   \code{"gaussian"} for continuous outcome(s) or \code{"binomial"} for binary
#'   outcome(s).
#' @param q number of outcome variables.
#' TODO change theta to include option for interactions
#' @param theta binary matrix with as many rows as predictors and as many
#'   columns as outcomes. A nonzero entry on row \eqn{i} and column \eqn{j}
#'   indicates that predictor \eqn{i} contributes to the definition of outcome
#'   \eqn{j}.
#' @param nu_xy vector of length \code{q} with expected proportion of predictors
#'   contributing to the definition of each of the \code{q} outcomes.
#' @param nu_int vector of length \code{q} with expected proportion of possible pairwise interactions
#'   contributing to the definition of each of the \code{q} outcomes.
#' @param hierarchy type of hierarchy between true predictors and considered pariwise interactions. 
#'   Possible values include "strong", only pairwise interactions between true predictor variables are considered,
#'   "weak", only pairwise interactions between one true predictor and one non-predictor variable are considered,
#'   "anti", only pairwise interaction between non-predictor variables are considered.
#' @param beta_abs vector defining the range of nonzero main effect regression coefficients
#'   in absolute values. If \code{continuous=FALSE}, \code{beta_abs} is the set
#'   of possible precision values. If \code{continuous=TRUE}, \code{beta_abs} is
#'   the range of possible precision values. Note that regression coefficients
#'   are re-scaled if \code{family="binomial"} to ensure that the desired
#'   concordance statistic can be achieved (see argument \code{ev_xy}) so they
#'   may not be in this range.
#' @param beta_int_abs vector defining the range of nonzero interaction effect regression coefficients
#'   in absolute values. If \code{beta_int_abs = NULL} \code{beta_abs} is used instead.
#'   If \code{continuous=FALSE}, \code{beta_int_abs} is the set of possible precision values.
#'   If \code{continuous=TRUE}, \code{beta_int_abs} is the range of possible precision values.
#'   Note that regression coefficients are re-scaled if \code{family="binomial"} to ensure that the desired
#'   concordance statistic can be achieved (see argument \code{ev_xy}) so they
#'   may not be in this range.
#' @param beta_sign vector of possible signs for regression coefficients.
#'   Possible inputs are: \code{1} for positive coefficients, \code{-1} for
#'   negative coefficients, or \code{c(-1, 1)} for both positive and negative
#'   coefficients.
#' @param continuous logical indicating whether to sample regression
#'   coefficients from a uniform distribution between the minimum and maximum
#'   values in \code{beta_abs} (if \code{continuous=TRUE}) or from proposed
#'   values in \code{beta_abs} (if \code{continuous=FALSE}).
#' @param ev_xy vector of length \code{q} with expected goodness of fit measures
#'   for each of the \code{q} outcomes. If \code{family="gaussian"}, the vector
#'   contains expected proportions of variance in each of the \code{q} outcomes
#'   that can be explained by the predictors. If \code{family="binomial"}, the
#'   vector contains expected concordance statistics (i.e. area under the ROC
#'   curve) given the true probabilities.

SimulateInteraction <- function(n = 100, pk = 10, xdata = NULL,
                                family = "gaussian", q = 1,
                                theta = NULL, nu_xy = 0.2,
                                nu_int = 1, hierarchy = "strong",
                                beta_abs = c(0.1, 1), beta_int_abs = NULL, beta_sign = c(-1, 1), continuous = TRUE,
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
  
  # TODO add option to implement own list of interactions and checks
  
  # Sampling pairwise interactions and interaction coefficients
  int <- SampleInteractions(xdata = xdata, theta = theta, q = q, nu_int = nu_int, hierarchy = hierarchy,
                            beta_abs = beta_abs, beta_int_abs = beta_int_abs, beta_sign = beta_sign, continuous = continuous)
  
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
  class(out) <- "simulation_regression"
  
  return(out)
}




SampleInteractions <- function(xdata, theta, q = NULL,
                               nu_int = 0.1,
                               beta_abs = c(0.1, 1), beta_int_abs = c(0.1, 1), beta_sign = c(-1, 1), continuous = TRUE, 
                               hierarchy = "strong"){
  
  # TODO remove for loop and replace with lapply
  
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
  
  ## instead of interaction table make into array ? potentially easier to sample beta
  
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
    
    # Defining interaction effect values
    
    if(!is.null(beta_int_abs)){
      beta_abs <- beta_int_abs
    }
    
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
    # ordering interactions (lower variable name first)
    
    int_table[[k]]$Var1 <- ordered(int_table[[k]]$Var1, colnames(xdata))
    int_table[[k]]$Var2 <- ordered(int_table[[k]]$Var2, colnames(xdata))
    
    # doesn't work?
    int_table[[k]] <- transform(int_table[[k]], Var1 = ifelse(Var1 > Var2, Var2, Var1), Var2 = ifelse(Var2 < Var1, Var1, Var2))
    int_table[[k]]$Var1 <- ordered(int_table[[k]]$Var1, 1:ncol(xdata))
    int_table[[k]]$Var2 <- ordered(int_table[[k]]$Var2, 1:ncol(xdata))
    
    levels(int_table[[k]]$Var1) <- colnames(xdata)
    levels(int_table[[k]]$Var2) <- colnames(xdata)
    
  }
  
  return(int_table) 
  
}