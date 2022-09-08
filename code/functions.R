



buster_cmax <- function(cmax, title){
  out <- log(cmax) - mean(log(cmax))
  
  barplot(height = out, ylab = "log(cmax) - mean(log(cmax))", col = "blue", main = title)
 
}



#linear mixed model with y = cmaxT
# x1 = subject nested in sequence
#x2 = treatment
#x3 = period
#x4 = sequence
#then plot cumulative confidence intervals from 5 subjects until end
#each iteration reports the last fitted value of the model and its confidence interval
mse <- function(sm){
  return(mean(resid(sm)^2))
}

#buster_conf <- function(data){
 # conflist <- list()
  #msevec <- rep(NA, nrow(data) - 9)
  #fitval <- rep(NA, nrow(data) - 9)
  #for(i in seq(10,nrow(data), by =2)){
  #  dat <- data[1:i,]
  #  mod <- lmer(log(cmax) ~ treat + sequence + period + (1|subject), data = dat)
  #  conflist[[i-9]] <- exp(confint(mod, c("treatT"), level = 0.9))
  #  msevec[i-9] <- mse(mod)
  #  fitval[i-9] <- exp(coef(mod)$subject$treatT[1])
  #  if(i == nrow(dat)){
  #    resi <- resid(mod)
  #  }
  #}
  #msevec <- na.omit(msevec)
  #fitval <- na.omit(fitval)
  #conflist <- conflist[-seq(2, length(conflist), by = 2)]
  #returnlist <- list(conflist, msevec, fitval, resi)
  #return(returnlist)
#}

buster_conf <- function(data){
  conflist <- list()
  msevec <- rep(NA, nrow(data) - 9)
  fitval <- rep(NA, nrow(data) - 9)
  for(i in seq(10,nrow(data), by =2)){
    dat <- data[1:i,]
    mod <- lm(log(cmax) ~ treat + sequence + period + subject, data = dat)
    conflist[[i-9]] <- exp(confint(mod, c("treatT"), level = 0.9))
    msevec[i-9] <- mse(mod)
    fitval[i-9] <- exp(coef(mod)[2])
    if(i == nrow(dat)){
      resi <- resid(mod)
    }
  }
  msevec <- na.omit(msevec)
  fitval <- na.omit(fitval)
  conflist <- conflist[-seq(2, length(conflist), by = 2)]
  returnlist <- list(conflist, msevec, fitval, resi)
  return(returnlist)
}





buster <- function(data){
  defaultW <- getOption("warn")
  options(warn=-1)
  {res <- buster_conf(data)
  
  devAskNewPage(ask = TRUE)
  buster_cmax(data$cmax[data$treat == "T"], title = "Buster Bar Plot Cmax Test")
  buster_cmax(data$cmax[data$treat == "R"], title = "Buster Bar Plot Cmax Reference")
  buster_confplot(res[[1]], res[[3]])
  buster_mseplot(res[[2]])
  buster_residplot(res[[4]][data$treat == "T"], title = "Residuals of Linear Model - Test")
  buster_residplot(res[[4]][data$treat == "R"], title = "Residuals of Linear Model - Reference")
  devAskNewPage(ask = FALSE)
  }
  options(warn = defaultW)
}







buster_confplot <- function(conflist, fitval){
  replwr <- rep(c(TRUE, FALSE), length(fitval)/2)
  repupr <- rep(c(FALSE, TRUE), length(fitval)/2)
  plot(0,0, xlim = c(5, length(fitval)+4), ylim = c(0.5,1.3), type = "n", xlab = "No. of Subjects in Analysis",
       ylab = "Confidence Interval", main = "Point estimates & 90%-CI Treat")
  lines(5:(length(fitval)+4), fitval, col = "blue", type = "pch", pch = 16)
  lines(5:(length(fitval)+4), unlist(conflist)[replwr], col = "red", type = "pch", pch = 16)
  lines(5:(length(fitval)+4), unlist(conflist)[repupr], col = "red", type = "pch", pch = 16)
  abline(h=0.8, col = "red")
  abline(h = 1, col= "blue")
  abline(h=1.25, col = "red")
  
}



buster_mseplot <- function(msevec){
  plot(x = 5:(length(msevec)+4), y = msevec, col = "purple", xlab = "No. of Subjects in Analysis",
       ylab = "MSE", pch = 16, main = "MSE of Linear Model")
}

buster_residplot <- function(resi, title){
  barplot(height = resi, col = "peachpuff2", xaxt = "n", ylab = "Residual(Cmax)", 
          main = title)
}

score32 <- function(v1, v2){
  #if(any(v1 == 0 & v1 == v2)){
   # v1[which(v1 == 0 & v1 == v2)] <- v1[which(v1 == 0 & v1 == v2)] + .Machine$double.xmin
    #v2[which(v2 == 0 & v1 == v2)] <- v2[which(v2 == 0 & v1 == v2)] + .Machine$double.xmin
  #}
  out <- 1/length(v1) * sum(abs(v1-v2)/(0.5 *((v1+v2)+.Machine$double.xmin)))
  
  return(out)
}


SaToWIB <- function(data){
  mat <- matrix(c(0,0,0,0),1,4)
  for(i in 1:(nrow(data)-1)){
    v1 <- sort(data[i, ], decreasing = TRUE)
    dilsum1<- sum(v1[1:3])
    data1 <- data[i,]/dilsum1
    for(j in (i+1):nrow(data)){
      v2 <- sort(data[j, ], decreasing = TRUE)
      dilsum2 <- sum(v2[1:3])
      data2 <- data[j,]/dilsum2
      
      mat <- rbind(mat, c(rownames(data)[i], rownames(data)[j], round(score32(data1, data2),5), 
                          round(dilsum1/dilsum2, 4)))
    }
  }
  return(mat[2:nrow(mat),])
}


buster_test <- function(data, alpha = 0.05){
 res <- buster_conf(data)
 point_est <- exp((res[[3]]+1)^2)
 #cummean <- rep(0, length(point_est))
 #for(i in 1:length(point_est)){
#   cummean[i] <- mean(point_est[max(c(1,i-12)):min(i,length(point_est)-7)])
# }
 #cummean <- cummean[(ceiling(length(point_est)/4)+1):length(point_est)]
 
 new <- point_est #- cummean
 
 kpss_test <- kpss.test(new)
 
 adf_test <- adf.test(new, alternative = "explosive")
 
 
 if(kpss_test$p.value < alpha & adf_test$p.value < 2 * alpha){
   out <- TRUE
 } else{
   out <- FALSE
 }
 return(out)
}


kpss_test <- function(data, alpha = 0.05){
  res <- buster_conf(data)
  point_est <- res[[3]]
  cummean <- rep(0, length(point_est))
  for(i in 1:length(point_est)){
     cummean[i] <- mean(point_est[max(c(1,i-7)):min(i,length(point_est)-7)])
   }
  #cummean <- cummean[(ceiling(length(point_est)/4)+1):length(point_est)]
  
  new <- point_est - cummean
  
  kpss_test <- kpss.test(new)
  
  
  
  if(kpss_test$p.value < alpha){
    out <- TRUE
  } else{
    out <- FALSE
  }
  return(out)
}

kpss_test <- function(data, alpha = 0.05){
  res <- buster_conf(data)
  point_est <- res[[3]] 
  cummean <- rep(0, length(point_est))
  for(i in 1:length(point_est)){
    cummean[i] <- mean(point_est[max(1, i-1):min(length(point_est), i+1)])
  }
  
  cummean[(length(point_est)*3/4):(length(point_est))] <- 0
  #cummean <- cummean[(ceiling(length(point_est)/4)+1):length(point_est)]
  
  new <- point_est - cummean
  
  kpss_test <- kpss.test(new)
  
  if(kpss_test$p.value < alpha){
    out <- TRUE
  } else{
    out <- FALSE
  }
  return(out)
}







cmaxdiff_test <- function(data, alpha = 0.05){
  point_est <- log(data$cmax[data$treat == "T"]) - mean(log(data$cmax[data$treat == "T"]))
  
  #point_est <- point_est[(ceiling(length(point_est)/3)+1):length(point_est)]
  
  test <- kpss.test(point_est)
  
  if(test$p.value < alpha){
    out <- TRUE
  } else{
    out <- FALSE
  }
  return(out)
}

library(strucchange)

strucc_test <- function(data, alpha = 0.05){
  res <- buster_conf(data)
  point_est <- res[[3]]

  ocus.res <- efp(point_est ~ 1, type = "OLS-CUSUM")
  test <- sctest(ocus.res)
  
  
  if(test$p.value < alpha){
    out <- TRUE
  } else{
    out <- FALSE
  }
  return(out)
}



confinttest <- function(data, lag = 3){
  res <- buster_conf(data)
  
  replwr <- rep(c(TRUE, FALSE), length(res[[3]])/2)
  repupr <- rep(c(FALSE, TRUE), length(res[[3]])/2)
  
  confints <- unlist(res[[1]])
  
  confdiff <- confints[repupr] - confints[replwr]
  
  diff_confdiff <- diff(confdiff)
  
  counter <- 0
  
  out <- FALSE
  for(i in  2:length(diff_confdiff)){
    
    if(diff_confdiff[i] < 0 & diff_confdiff[i-1] > 0){
      counter <- 0
    }
    if(diff_confdiff[i] > 0){
      counter <- counter + 1
    }
    
    if(counter >= lag){
      out <- TRUE
      break 
    } 
  }

  return(out)
  
}

confwidth_test <- function(data, alpha = 0.05){
  res <- buster_conf(data)
  replwr <- rep(c(TRUE, FALSE), length(res[[3]])/2)
  repupr <- rep(c(FALSE, TRUE), length(res[[3]])/2)
  
  confints <- unlist(res[[1]])
  
  confdiff <- confints[repupr] - confints[replwr]
  
  test <- kpss.test(confdiff)
  
  if(test$p.value < alpha){
    out <- TRUE
  } else{
    out <- FALSE
  }
  return(out)
}


confwidthstrucc_test <- function(data, alpha = 0.05){
  res <- buster_conf(data)
  replwr <- rep(c(TRUE, FALSE), length(res[[3]])/2)
  repupr <- rep(c(FALSE, TRUE), length(res[[3]])/2)
  
  confints <- unlist(res[[1]])
  
  confdiff <- confints[repupr] - confints[replwr]
  
  ocus.res <- efp(confdiff ~ 1, type = "OLS-CUSUM")
  test <- sctest(ocus.res)
  
  
  if(test$p.value < alpha){
    out <- TRUE
  } else{
    out <- FALSE
  }
  return(out)
}



adfexp_test <- function(data, alpha = 0.05){
  res <- buster_conf(data)
  point_est <- exp((res[[3]]+1)^4)
  point_est <- point_est[(ceiling(length(point_est)/3)+1):length(point_est)]
  
  test <- adf.test(point_est, alternative = "explosive")
  
  if(test$p.value < alpha){
    out <- TRUE
  } else{
    out <- FALSE
  }
  return(out)
}


adf_test <- function(data, alpha = 0.05){
  res <- buster_conf(data)
  point_est <- res[[3]]
  point_est <- point_est[(ceiling(length(point_est)/3)+1):length(point_est)]
  
  test <- adf.test(point_est, alternative = "stationary")
  
  if(test$p.value > alpha){
    out <- TRUE
  } else{
    out <- FALSE
  }
  return(out)
}





