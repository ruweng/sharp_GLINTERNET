#### GLINTERNET Performance and Visuals ####

### Heatmap colour gradient ----

colgradient <- colorRampPalette(c("#6bcbb8", "#ffc500"))

### Selection proportion Performance ----

plot.selprop_performance <- function(x,
                                       col = c("#6bcbb8","#ff0000", "#ffc500","#ff5c02"),
                                       col.axis = NULL,
                                       col.thr = "#A0AAA3",
                                       lty.thr = 2,
                                       n_predictors = NULL,
                                       theta = theta,
                                       ...) {
  
  
  # Storing extra arguments
  extra_args <- list(...)
  
  # Defining default settings
  if (!"type" %in% names(extra_args)) {
    extra_args$type <- "h"
  }
  if (!"xlab" %in% names(extra_args)) {
    extra_args$xlab <- ""
  }
  if (!"ylab" %in% names(extra_args)) {
    extra_args$ylab <- "Selection proportion"
  }
  if (!"las" %in% names(extra_args)) {
    extra_args$las <- 2
  }
  if (!"ylim" %in% names(extra_args)) {
    extra_args$ylim <- c(0, 1)
  }
  if (!"cex.lab" %in% names(extra_args)) {
    extra_args$cex.lab <- 1.2
  }
  if (!"cex.axis" %in% names(extra_args)) {
    extra_args$cex.axis <- 1
  }
  
  # Checking inputs
  if (length(col) != 4) {
    col <- rep(col[1], )
  }
  if (is.null(col.axis)) {
    col.axis <- col
  } else {
    if (length(col.axis) != 4) {
      col.axis <- rep(col.axis[1], 4)
    }
  }
  
  
  # Extracting selection proportions
  selprop <- SelectionProportions(x)
  
  # Defining colours
  mycolours <- ifelse(SelectedVariables(x) == 1, 
                      yes = ifelse(SelectedVariables(x) == theta, yes = col[1], no = col[2]),
                      no = ifelse(SelectedVariables(x) == theta, yes = col[3], no = col[4]))
  mycolours_axis <- mycolours
  
  # Re-ordering by decreasing selection proportion
  ids <- sort.list(selprop, decreasing = TRUE)
  if (is.null(n_predictors)) {
    n_predictors <- length(ids)
  }
  n_predictors <- min(n_predictors, length(ids))
  ids <- ids[1:n_predictors]
  selprop <- selprop[ids]
  mycolours <- mycolours[ids]
  mycolours_axis <- mycolours_axis[ids]
  
  # Extracting relevant extra arguments
  tmp_extra_args <- extra_args
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("xaxt", "col")]
  
  # Making plot
  do.call(base::plot, args = c(
    list(
      x = selprop,
      xaxt = "n",
      col = mycolours
    ),
    tmp_extra_args
  ))
  
  # Extracting relevant extra arguments
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = graphics::axis)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("side", "at", "labels", "col.axis", "col.ticks", "type")]
  
  # Adding axis
  for (i in 1:length(selprop)) {
    do.call(graphics::axis, args = c(
      list(
        side = 1,
        at = i,
        labels = names(selprop)[i],
        col.axis = mycolours_axis[i],
        col.ticks = mycolours_axis[i]
      ),
      tmp_extra_args
    ))
  }
  graphics::abline(h = Argmax(x)[1, 2], lty = lty.thr, col = col.thr)
  
  # Extracting performance measures
  perf <- SelectionPerformance(x, theta)
  
  # Adding legend 
  
  legend(x = "topright",
         lty = c(4,6),
         col= append(col,col.thr), 
         box.lty= 0,
         cex = 0.6,
         legend=c(paste0("True Positives N(",perf$TP,")"), 
                  paste0("False Positives N(",perf$FP,")"),
                  paste0("True Negatives N(",perf$TN,")"),
                  paste0("False Negatives N(",perf$FN,")"),
                  expression(hat(pi))))
}

### Stability Paths ----

plot.stab_paths <- function(stab_int, sim_int){
  
  # custom functions
  addTrans <- function(color,trans)
  {
    # This function adds transparancy to a color.
    # Define transparancy with an integer between 0 and 255
    # 0 being fully transparant and 255 being fully visable
    # Works with either color and trans a vector of equal length,
    # or one of the two of length 1.
    
    if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
    if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
    if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
    
    num2hex <- function(x)
    {
      hex <- unlist(strsplit("0123456789ABCDEF",split=""))
      return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
    }
    rgb <- rbind(col2rgb(color),trans)
    res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
    return(res)
  }
  
  # create dataframe of stability paths
  
  Pi <- stab_int$P
  
  Lambda <- stab_int$Lambda
  
  index_lambda <- 1:length(stab_int$Lambda)
  
  selprop <- t(sapply(index_lambda, FUN = function(x){SelectionProportions(stab_int, x)}))
  
  theta <- SelectedVariables(stab_int)
  
  theta_star <- sim_int$theta
  
  beta_true <- sim_int$beta
  
  stab_paths <- as.data.frame(cbind(Lambda, selprop))
  
  # define margins and layout
  
  par(oma = c(2,2,2,2))
  par(mar=c(4,4,4,10), mfrow = c(1,1), las = 1)
  
  # empty plot 
  
  plot(NA,
       xlim = c(min(log(Lambda)), ceiling(max(log(Lambda), na.rm = T))),
       ylim = c(0, 1),
       frame.plot =  FALSE, 
       axes = T,
       xlab = NA,
       ylab = NA,
       xaxt = "n",
       yaxt = "n",
       xaxs = "i",
       yaxs = "i")
  
  # add axes
  
  axis(side = 1, at = seq(floor(round(min(log(Lambda)), digits = 2)), ceiling(max(log(Lambda))), by = 1))
  axis(side = 2, at = seq(0,1, by = 0.1), las=1)
  
  # Add titles
  
  title(
    adj = 0,
    line = 2,
    cex.axis=0.5, 
    cex.main=1, 
    main = expression(paste("Selection proportion as a Function of ", lambda)))
  
  title( 
    cex.lab=0.9, 
    ylab = "Selection proportion")
  
  title(
    cex.lab=0.9,
    adj = 1,
    xlab = expression(paste("log(",lambda,")")))
  
  # add abline for calibrated parameters 
  abline(v = log(Argmax(stab_int)[1]), lty = 2, col = addTrans("#A0AAA3", 255))
  abline(h = Argmax(stab_int)[2], lty = 2, col = addTrans("#A0AAA3", 255))
  lines(log(Lambda), Pi, col = addTrans("#A0AAA3", 150))
  
  # add stability paths 
  for (i in 2:ncol(stab_paths)){
    if (theta[colnames(stab_paths)[i]] == 1 & theta[colnames(stab_paths)[i]] == theta_star[colnames(stab_paths)[i],]){
      lines(log(Lambda), stab_paths[,i], col = addTrans("#6bcbb8", trans =  255), lwd = beta_true[colnames(stab_paths)[i],] + 1)
      text(x = log(Lambda[max(which(unique(stab_paths[,i]) != 1))]), y = stab_paths[,i][max(which(unique(stab_paths[,i]) != 1))], labels =  colnames(stab_paths)[i],
           cex=0.65, col = addTrans("#6bcbb8", trans =  255))
    } else if(theta[colnames(stab_paths)[i]] == 0 & theta[colnames(stab_paths)[i]] == theta_star[colnames(stab_paths)[i],]){
      lines(log(Lambda), stab_paths[,i], col = addTrans("#ffc500", trans = 70), lty = 2, lwd = beta_true[colnames(stab_paths)[i],] + 1)
    } else if(theta[colnames(stab_paths)[i]] == 1 & theta[colnames(stab_paths)[i]] != theta_star[colnames(stab_paths)[i],]){
      lines(log(Lambda), stab_paths[,i], col = addTrans("#ff0000", trans = 255))
      text(x = log(Lambda[max(which(unique(stab_paths[,i]) != 1))]), y = stab_paths[,i][max(which(unique(stab_paths[,i]) != 1))], labels =  colnames(stab_paths)[i],
           cex=0.65, col = addTrans("#ff0000", trans =  255))
    } else{ 
      lines(log(Lambda), stab_paths[,i], col = addTrans("#ff5c02", trans = 255), lty = 2, lwd = beta_true[colnames(stab_paths)[i],] + 1)
      text(x = log(Lambda[max(which(unique(stab_paths[,i]) != 1))]), y = stab_paths[,i][max(which(unique(stab_paths[,i]) != 1))], labels =  colnames(stab_paths)[i],
           cex=0.65, col = addTrans("#ff5c02", trans =  255))
    }
  }
  
  # add gridlines 
  
  grid(ny = 10)
  
  # add stability as background?
  
  legend(x = "topright",
         legend = c("True Positives", "False Positives", "True Negatives", "False Negatives", substitute(hat(pi)), substitute(paste(hat(lambda), " , ", hat(pi)))),
         col = c("#6bcbb8","#ff0000", "#ffc500","#ff5c02", "#A0AAA3", "#A0AAA3"),
         lty=c(1,1,2,2,1,2),
         cex=1,
         inset = c(-0.2, 0),
         horiz= F,
         xpd = T,
         bty = "n",
         seg.len= 1,
         x.intersp = 0.5,
         y.intersp = 1.5)
  
}

### Grid Performance ----

grid_performance <- function(stab_int, sim_int){
  
  # define grid of parameter indices
  
  index_lambda = 1:nrow(stab_int$S_2d) 
  index_pi = 1:ncol(stab_int$S_2d)
  
  index <- expand_grid(index_lambda, index_pi)
  
  # extract stably selected variables for each parameter combination
  
  theta_grid <- apply(index, 1, FUN = function(x){Stable(stab_int, x)})
  
  # extract stability score for each parameter combination
  
  stability_score <- as.vector(t(stab_int$S_2d))
  
  # extract parameter combinations
  
  params_grid <- expand_grid(stab_int$Lambda,stab_int$params$pi_list)
  colnames(params_grid) <- c("lambda","pi")
  
  # calculate selection performance per parameter combination
  
  performance <- as.data.frame(apply(mapply("[",(apply(theta_grid, 2, FUN = function(x){SelectionPerformance(theta = x, theta_star = sim_int$theta)}))),1, unlist))
  
  str(performance)
  
  # combine parameter combination and performance 
  
  grid_performance <- cbind(params_grid, stability_score, performance)
  
  # add logical of max stability score
  
  grid_performance$max_stability <- as.factor(ifelse(grid_performance$stability_score == max(grid_performance$stability_score,  na.rm = TRUE),1,0))
  
  return(grid_performance)
}

plot.grid_performance <- function(grid_performance){
  
  # custom functions
  addTrans <- function(color,trans)
  {
    # This function adds transparancy to a color.
    # Define transparancy with an integer between 0 and 255
    # 0 being fully transparant and 255 being fully visable
    # Works with either color and trans a vector of equal length,
    # or one of the two of length 1.
    
    if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
    if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
    if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
    
    num2hex <- function(x)
    {
      hex <- unlist(strsplit("0123456789ABCDEF",split=""))
      return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
    }
    rgb <- rbind(col2rgb(color),trans)
    res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
    return(res)
  }
  
  # define margins and layout
  
  par(oma=c(2,2,2,2))
  par(mar=c(4,4,4,4),
      mfrow=c(1,2))
  
  # F1 vs Stability ----
  
  # create empty plot
  
  plot(
    NA,
    xlim = c(0, max(grid_performance$stability_score, na.rm = T)),
    ylim = c(0, 1),
    frame.plot =  FALSE, 
    axes = T,
    xlab = NA,
    ylab = NA,
    xaxt = "n",
    yaxt = "n",
    xaxs = "i",
    yaxs = "i"
  )
  
  # add axes
  
  axis(side = 1, at = append(seq(0,max(grid_performance$stability_score, na.rm = T), by = 10),ceiling(max(grid_performance$stability_score, na.rm = T))))
  axis(side = 2, at = seq(0,1, by = 0.1), las=1)
  
  # Add titles
  
  title(
    adj = 0,
    line = 2,
    cex.axis=0.5, 
    cex.main=1, 
    main = expression("A) as a function of the stability score"))
  
  title( 
    cex.lab=0.9, 
    ylab = "F1 score")
  
  title(
    cex.lab=0.9,
    adj = 1,
    xlab = "Stability score",)
  
  # add grid lines
  
  grid(ny = 10)
  
  
  # add points
  
  points( x = grid_performance$stability_score[grid_performance$max_stability == 0],
          y = grid_performance$F1_score[grid_performance$max_stability == 0],
          cex = 0.5,
          pch = 17,
          col = addTrans("#ffc500", trans= 70))
  abline(v = grid_performance$stability_score[grid_performance$max_stability == 1], 
         lty = 2, 
         col = addTrans("#A0AAA3", trans= 255)) 
  abline(h = grid_performance$F1_score[grid_performance$max_stability == 1], 
         lty = 2,
         col = addTrans("#A0AAA3", trans= 255)) 
  points( x = grid_performance$stability_score[grid_performance$max_stability == 1],
          y = grid_performance$F1_score[grid_performance$max_stability == 1],
          cex = 1,
          pch = 17,
          xpd= T, 
          col = c("#6bcbb8"))
  
  # add legend
  
  # F1 vs Lambda ----
  
  # create empty plot
  
  plot(
    NA,
    xlim = c(min(log(grid_performance$lambda)), ceiling(max(log(grid_performance$lambda), na.rm = T))),
    ylim = c(0, 1),
    frame.plot =  FALSE, 
    axes = T,
    xlab = NA,
    ylab = NA,
    xaxt = "n",
    yaxt = "n",
    xaxs = "i",
    yaxs = "i"
  )
  
  # add axes
  
  axis(side = 1, at = seq(floor(round(min(log(grid_performance$lambda)), digits = 2)), ceiling(max(log(grid_performance$lambda))), by = 1))
  axis(side = 2, at = seq(0,1, by = 0.1), las=1)
  
  # Add titles
  
  title(
    adj = 0,
    line = 2,
    cex.axis=0.5, 
    cex.main=1, 
    main = expression(paste("B) as a function of ", lambda)))
  
  title( 
    cex.lab=0.9, 
    ylab = "F1 score")
  
  title(
    cex.lab=0.9,
    adj = 1,
    xlab = expression(paste("log(",lambda,")")))
  
  # add grid lines
  
  grid(ny = 10)
  
  
  # add points
  
  points( x = log(grid_performance$lambda[grid_performance$max_stability == 0]),
          y = grid_performance$F1_score[grid_performance$max_stability == 0],
          cex = 0.5,
          pch = 17,
          col = addTrans("#ffc500", trans= 70))
  abline(v = log(grid_performance$lambda[grid_performance$max_stability == 1]), 
         lty = 2, 
         col = addTrans("#A0AAA3", trans= 255)) 
  abline(h = grid_performance$F1_score[grid_performance$max_stability == 1], 
         lty = 2,
         col = addTrans("#A0AAA3", trans= 255)) 
  points( x = log(grid_performance$lambda[grid_performance$max_stability == 1]),
          y = grid_performance$F1_score[grid_performance$max_stability == 1],
          cex = 1,
          pch = 17,
          col = c("#6bcbb8"))
  
  
  # add title 
  
  mtext(text = substitute(paste("Model performance across parameter combinations (", lambda, "," , pi, ")")),
        adj = 0, 
        side = 3,
        line = 0,
        outer = TRUE, 
        cex = 1.2)
  
  # add legend
  
}

### Cross-validation CalibrationPlot ----

plot.cv_calibration <- function(cv, model = "glinternet"){
  
  # custom functions
  addTrans <- function(color,trans)
  {
    # This function adds transparancy to a color.
    # Define transparancy with an integer between 0 and 255
    # 0 being fully transparant and 255 being fully visable
    # Works with either color and trans a vector of equal length,
    # or one of the two of length 1.
    
    if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
    if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
    if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
    
    num2hex <- function(x)
    {
      hex <- unlist(strsplit("0123456789ABCDEF",split=""))
      return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
    }
    rgb <- rbind(col2rgb(color),trans)
    res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
    return(res)
  }
  
  # define variables
  
  if(model == "glinternet"){
    Lambda = cv$lambda
    cvErr = cv$cvErr
    cvErrStd = cv$cvErrStd
    cvErrmax = cvErr + cvErrStd
    cvErrmin = cvErr - cvErrStd
    lambdahat = cv$lambdaHat
    lambdahat1Std = cv$lambdaHat1Std
  } else if(model == "lasso"){
    Lambda = cv$lambda
    cvErr = cv$cvm
    cvErrStd = cv$cvsd
    cvErrmax = cv$cvup
    cvErrmin = cv$cvlo
    lambdahat = cv$lambda.min
    lambdahat1Std = cv$lambda.1se
  }
  
  
  # define margins and layout
  
  par(oma = c(1,1,1,1))
  par(mar=c(4,4,4,8), mfrow = c(1,1), las = 1)
  
  # empty plot 
  
  plot(NA,
       xlim = c(min(log(Lambda), na.rm = T), max(log(Lambda), na.rm = T)),
       ylim = c(0, max(cvErrmax)),
       frame.plot =  FALSE, 
       axes = T,
       xlab = NA,
       ylab = NA,
       xaxt = "n",
       yaxt = "n",
       xaxs = "i",
       yaxs = "i")
  
  # add axes
  
  axis(side = 1, at = seq(floor(round(min(log(Lambda)), digits = 2)), ceiling(max(log(Lambda))), by = 1))
  axis(side = 2, at = seq(0,max(cvErrmax), by = 0.5))
  
  # Add titles
  
  title(
    adj = 0,
    line = 2,
    cex.axis=0.5, 
    cex.main=1, 
    main = expression(paste("Cross-validation Error as a Function of ", lambda)))
  
  title( 
    cex.lab=0.9, 
    ylab = "Cross-validation Error")
  
  title(
    cex.lab=0.9,
    adj = 1,
    xlab = expression(paste("log(",lambda,")")))
  
  # add abline for calibrated parameters
  abline(v = log(lambdahat), lty = 2, col = addTrans("#ff0000", 255))
  abline(v = log(lambdahat1Std), lty = 2, col = addTrans("#ffc500", 255))
  
  # plot cross-validation error
  lines(log(Lambda), cvErr)
  arrows(x0=log(Lambda), y0=cvErrmin, x1=log(Lambda), y1=cvErrmax, code=3, angle=90, col = addTrans("#A0AAA3", 150), length=0.05)
  
  # add legend

  legend(x = "topright",
         legend = c(substitute(paste(hat(lambda))), expression(paste(hat(lambda), " 1",sigma))),
         col = c("#ff0000", "#ffc500"),
         lty=c(2,2),
         cex=0.8,
         inset = c(-0.2, 0),
         horiz= F,
         xpd = T,
         bty = "n",
         seg.len= 1,
         x.intersp = 0.5,
         y.intersp = 1.5)
}

### Regularization Paths ----

plot.reg_paths <- function(reg_beta = NULL,
                           cv = cv,
                           calibration = "lambdahat_sd",
                           sim_int = sim_int){
  
  # custom functions
  addTrans <- function(color,trans)
  {
    # This function adds transparancy to a color.
    # Define transparancy with an integer between 0 and 255
    # 0 being fully transparant and 255 being fully visable
    # Works with either color and trans a vector of equal length,
    # or one of the two of length 1.
    
    if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
    if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
    if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
    
    num2hex <- function(x)
    {
      hex <- unlist(strsplit("0123456789ABCDEF",split=""))
      return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
    }
    rgb <- rbind(col2rgb(color),trans)
    res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
    return(res)
  }
  
  # retrieve plotting items
  
  if (is.null(reg_beta)){
    
    if(calibration == "lambdahat"){
      
      coef_lambdahat <- coef(cv, s = "lambda.min")
      
      theta <- ifelse(as.matrix(coef_lambdahat) == 0, 0, 1)[-1,]
      
      lambda_calibrate <- cv$lambda.min
      
    } else if (calibration == "lambdahat_sd"){
      
      coef_lambdahat_sd <- coef(cv, s = "lambda.1se")
      
      theta <- ifelse(as.matrix(coef_lambdahat_sd) == 0, 0, 1)[-1,]
      
      lambda_calibrate <- cv$lambda.1se
      
    }
    
    Lambda <- cv$lambda
    
    beta <- t(as.matrix(coef(cv, s = cv$lambda)))[,-1]
    
  } else {
  
  numLevels <- apply(sim_int$xdata,2, function(x) {
    if (is.factor(x))
      length(levels(x))
    else {1}})
  
  if(calibration == "lambdahat"){
    
    lambda_calibrate <- cv_int$lambdaHat
    
    glinternet_lambdahat <- glinternet::glinternet(sim_int$xdata, sim_int$ydata, lambda = cv_int$lambdaHat, numLevels = numLevels)
    
    theta <- glinternet.theta(mymodel = glinternet_lambdahat, xdata = sim_int$xdata)
    
  } else if (calibration == "lambdahat_sd"){
    
    lambda_calibrate <- cv_int$lambdaHat1Std
    
    glinternet_lambdahat_sd<- glinternet::glinternet(sim_int$xdata, sim_int$ydata, lambda = cv_int$lambdaHat1Std, numLevels = numLevels)
    
    theta <- glinternet.theta(mymodel = glinternet_lambdahat_sd, xdata = sim_int$xdata)
    
  }
  
  Lambda <- cv_int$lambda
  
  beta <- reg_beta
  
  }
  
  theta_star <- sim_int$theta
  
  beta_true <- sim_int$beta
  
  # define margins and layout
  
  par(oma = c(2,2,2,2))
  par(mar=c(4,4,4,10), mfrow = c(1,1), las = 1)
  
  # empty plot 
  
  plot(NA,
       xlim = c(min(log(Lambda)), ceiling(max(log(Lambda), na.rm = T))),
       ylim = c(floor((min(beta))),ceiling(max(beta))),
       frame.plot =  FALSE, 
       axes = T,
       xlab = NA,
       ylab = NA,
       xaxt = "n",
       yaxt = "n",
       xaxs = "i",
       yaxs = "i")
  
  # add axes
  
  axis(side = 1, at = seq(floor(round(min(log(Lambda)), digits = 2)), ceiling(max(log(Lambda))), by = 1))
  axis(side = 2, at = seq(floor((min(beta))),ceiling(max(beta)), by = 0.2), las=1)
  
  # Add titles
  
  title(
    adj = 0,
    line = 2,
    cex.axis=0.5, 
    cex.main=1, 
    main = expression(paste("Coefficients as a Function of ", lambda)))
  
  title( 
    cex.lab=0.9, 
    ylab = "Coefficients")
  
  title(
    cex.lab=0.9,
    adj = 1,
    xlab = expression(paste("log(",lambda,")")))
  
  # add abline for calibrated parameters 
  abline(v = log(lambda_calibrate), lty = 2, col = addTrans("#A0AAA3", 255))
  abline(h = 0, lty = 1, col = addTrans("black", 255))
  
  # add stability paths 
  for (i in 1:ncol(beta)){
    if (theta[colnames(beta)[i]] == 1 & theta[colnames(beta)[i]] == theta_star[colnames(beta)[i],]){
      lines(log(Lambda), beta[,i], col = addTrans("#6bcbb8", trans =  255), lwd = beta_true[colnames(beta)[i],] + 1)
      text(x = log(Lambda[max(which(unique(beta[,i]) != 1))]) + 0.5, y = beta[,i][max(which(unique(beta[,i]) != 1))], labels =  colnames(beta)[i],
           cex=0.65, col = addTrans("#6bcbb8", trans =  255))
      } else if(theta[colnames(beta)[i]] == 0 & theta[colnames(beta)[i]] == theta_star[colnames(beta)[i],]){
      lines(log(Lambda), beta[,i], col = addTrans("#ffc500", trans = 150), lty = 2, lwd = beta_true[colnames(beta)[i],] + 1)
      } else if(theta[colnames(beta)[i]] == 1 & theta[colnames(beta)[i]] != theta_star[colnames(beta)[i],]){
      lines(log(Lambda), beta[,i], col = addTrans("#ff0000", trans = 150))
      } else{ 
      lines(log(Lambda), beta[,i], col = addTrans("#ff5c02", trans = 255), lty = 2, lwd = beta_true[colnames(beta)[i],] + 1)
      text(x = log(Lambda[max(which(unique(beta[,i]) != 1))])  + 0.5, y = beta[,i][max(which(unique(beta[,i]) != 1))], labels =  colnames(beta)[i],
           cex=0.65, col = addTrans("#ff5c02", trans =  255))
    }
  }
  
  # add gridlines 
  
  grid(ny = 10)
  
  # add legend
  
  if (calibration == "lambdahat"){
    calibration_legend <- substitute(paste(hat(lambda)))
  } else if (calibration == "lambdahat_sd"){
    calibration_legend <- expression(paste(hat(lambda), " 1",sigma))
  }

  legend(x = "topright",
         legend = c("True Positives", "False Positives", "True Negatives", "False Negatives", calibration_legend),
         col = c("#6bcbb8","#ff0000", "#ffc500","#ff5c02", "#A0AAA3"),
         lty=c(1,1,2,2,2),
         cex=1,
         inset = c(-0.2, 0),
         horiz= F,
         xpd = T,
         bty = "n",
         seg.len= 1,
         x.intersp = 0.5,
         y.intersp = 1.5)
  
}
  


### Combine Performance Scripts ####

combine_performance <- function(){
  performance <- list.files(pattern = ".rds") %>%
    map(readRDS) %>% 
    bind_rows()
}

### Plot asymptotic performance ####

plot.asymptotic_performance <- function(performance, param = "ev_xy"){
  
  # custom functions
  addTrans <- function(color,trans)
  {
    # This function adds transparency to a color.
    # Define transparency with an integer between 0 and 255
    # 0 being fully transparent and 255 being fully visible
    # Works with either color and trans a vector of equal length,
    # or one of the two of length 1.
    
    if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
    if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
    if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
    
    num2hex <- function(x)
    {
      hex <- unlist(strsplit("0123456789ABCDEF",split=""))
      return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
    }
    rgb <- rbind(col2rgb(color),trans)
    res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
    return(res)
  }
  
  # define margins and layout
  par(oma=c(2,2,0,0))
  par(mar=c(5,5,5,5),
      mfrow=c(2,2))
  
  # define colors
  color = c("#ff0000", "#ffc500","#6bcbb8","#ff5c02")
  
  # plot for each number of variables pk
  pk <- as.factor(unlist(performance$pk))
  
  for (i in levels(pk)){
    
    # split dataset by pk
    performance_split <- performance[which(pk == i),]
    
    # retrieve parameters
    param_tmp <- unlist(performance_split[param])
    F1_tmp <- unlist(performance_split["F1_score"])
    recall_tmp <- unlist(performance_split["recall"])
    precision_tmp <- unlist(performance_split["precision"])
    model_tmp <- unlist(performance_split["model"])
    
    performance_tmp <- as.data.frame(cbind(model_tmp, param_tmp, F1_tmp, recall_tmp, precision_tmp)) %>% 
      mutate(model_tmp = as.factor(model_tmp),
             param_tmp = as.factor(param_tmp),
             F1_tmp = as.numeric(F1_tmp),
             recall_tmp = as.numeric(recall_tmp),
             precision_tmp = as.numeric(precision_tmp)
      )
    
    ## plot F1 as a function of param
    # create plot
    plot(NA,
         xlim = c(min(as.numeric(levels(performance_tmp$param_tmp))), max(as.numeric(levels(performance_tmp$param_tmp)))),
         ylim = c(0, max(performance_tmp$F1_tmp)),
         frame.plot =  FALSE, 
         axes = T,
         xlab = NA,
         ylab = NA,
         xaxt = "n",
         yaxt = "n",
         xaxs = "i",
         yaxs = "i")
    
    # add axes
    axis(side = 1, at = seq(0, max(as.numeric(levels(performance_tmp$param_tmp))), by = abs(as.numeric(levels(performance_tmp$param_tmp)[1])-as.numeric(levels(performance_tmp$param_tmp)[2]))))
    axis(side = 2, at = seq(0, max(performance_tmp$F1_tmp), by = 0.1), las=1)
    
    # create titles
    title(
      adj = 0,
      line = 2,
      cex.main=1, 
      main = paste0("F1 score as a function of ", param, " for pk = ", i))
    title( 
      cex.lab=1, 
      ylab = "F1 score")
    title(
      cex.lab=1,
      adj = 1,
      xlab = param)
    
    # compute summary statistics
    
    quibble <- function(x, q = c(0.25, 0.5, 0.75)) {
      tibble(x = quantile(x, q), q = q)
    }
    performance_summary <- performance_tmp %>% 
      group_by(model_tmp, param_tmp) %>% 
      reframe(F1 = quantile(F1_tmp, c(0.25, 0.5, 0.75)),
              Recall = quantile(recall_tmp, c(0.25, 0.5, 0.75)),
              Precision = quantile(precision_tmp, c(0.25, 0.5, 0.75)),
              q = c(0.25, 0.5, 0.75))
    
    # plot median lines and IQR ribbons
    for (k in 1:4){
      level <- levels(performance_summary$model_tmp)[k]
      F1_median <- performance_summary[which(performance_summary$model_tmp == level & performance_summary$q == 0.5), "F1"]
      F1_q1 <- performance_summary[which(performance_summary$model_tmp == level & performance_summary$q == 0.25), "F1"]
      F1_q3 <- performance_summary[which(performance_summary$model_tmp == level & performance_summary$q == 0.75), "F1"]
      param_summary <- performance_summary[which(performance_summary$model_tmp == level & performance_summary$q == 0.5), "param_tmp"]
      graphics::lines(as.numeric(unlist(F1_median)) ~ as.numeric(levels(unlist(param_summary))), lwd=1, col = color[[k]])
      graphics::polygon(c(as.numeric(levels(unlist(param_summary))), as.numeric(rev(levels(unlist(param_summary))))), c(as.numeric(unlist(F1_q1)),as.numeric(rev(unlist(F1_q3)))), col = addTrans(color[[k]], 20), border=NA)
    }
    
    ## plot average F1 per model
    # create plot
    plot(NA,
         xlim = c(0,length(levels(performance_tmp$model_tmp))),
         ylim = c(0, max(performance_tmp$F1_tmp)),
         frame.plot =  FALSE, 
         axes = T,
         xlab = NA,
         ylab = NA,
         xaxt = "n",
         yaxt = "n",
         xaxs = "i",
         yaxs = "i")
    
    # add axes
    axis(side = 1, at = seq(0.2,length(levels(performance_tmp$model_tmp))-0.8), labels = c("Cross-validation", "Cross-validation SD", "Stability Z-score", "Stability NLL"), las = 2, cex.axis = 0.7)
    axis(side = 2, at = seq(0, max(performance_tmp$F1_tmp), by = 0.1), las=1, cex.axis = 0.7)
    
    # create titles
    title(
      adj = 0,
      line = 2,
      cex.main=1, 
      main = paste0("Average F1 score across simulations for pk = ", i))
    title( 
      cex.lab=1, 
      ylab = "F1 score")
    
    
    
    # create summary statistics
    
    performance_mean <- performance_tmp %>% 
      group_by(model_tmp) %>% 
      summarise(F1_mean = mean(F1_tmp), 
                F1_sd = sd(F1_tmp))
    
    # plot means with standard deviation error bars
    
    points(seq(0.2,length(levels(performance_tmp$model_tmp))-0.8), performance_mean$F1_mean, col = color)
    segments(seq(0.2,length(levels(performance_tmp$model_tmp))-0.8), performance_mean$F1_mean - performance_mean$F1_sd, seq(0.2,length(levels(performance_tmp$model_tmp))-0.8), performance_mean$F1_mean + performance_mean$F1_sd, col = color)
    
  }
  
}