sim_data <- function(n, mean, sd, treat_effect, seed, manipulation = TRUE){

  set.seed(seed)
  logcmaxR <- exp(mean + rnorm(n, sd = sd))

  set.seed(seed+1)
  logcmaxT <- exp(mean + treat_effect  + rnorm(n, sd = sd))

  sortR <- sort(logcmaxR)

  sortT <- sort(logcmaxT)

  simdata_rt <- cbind(sortR, sortT)

  simdatadf <- data.frame(simdata_rt)
  colnames(simdatadf) <- c("R", "T")

  set.seed(seed+2)
  simdatadf <- simdatadf[sample(nrow(simdatadf)),]
  
  if(manipulation){
    simdatadf$ratios <- simdatadf$T / simdatadf$R
    sortdf <- simdatadf[1:ceiling(nrow(simdatadf)*(2/3)),]
    sorteddf <- sortdf[order(sortdf$ratios),]
    
    if(treat_effect < 0){
      candidates <- sorteddf[1:ceiling(nrow(sorteddf)/2),]
      candidates[1:ceiling(nrow(candidates)/4), ]$T <- candidates[1:ceiling(nrow(candidates)/4), ]$T * 0.5
      
    } else{
      candidates <- sorteddf[ceiling(nrow(sorteddf)/2+1): nrow(sorteddf), ]
      candidates[1:(ceiling(nrow(candidates)/4)), ]$T <- candidates[1:(ceiling(nrow(candidates)/4)), ]$T * (1/0.5)
    }
    tmp <- candidates$T
    candidates$T <- candidates$R
    candidates$R <- tmp

    
    simdatadf[(nrow(simdatadf)-nrow(candidates)+1):nrow(simdatadf), ] <- candidates
  }
  
  

  set.seed(34)
  seq_rand <- sample(c("RT","TR"),n,replace=T)

  sim_test <- simdatadf$T
  sim_ref <- simdatadf$R

  sub_R <- cbind(as.character(1:n),as.numeric(sim_ref),"R",seq_rand)
  colnames(sub_R) <- c("subject", "cmax","treat","sequence")

  sub_T <- cbind(as.character(1:n),as.numeric(sim_test),"T",seq_rand)
  colnames(sub_R) <- c("subject", "cmax","treat","sequence")

  simdf <- rbind(sub_R,sub_T)
  simdf <- data.frame(simdf)
  simdf$period <- 0
  simdf[which(simdf$sequence=="RT" & simdf$treat=="R"),]$period <- "1"
  simdf[which(simdf$sequence=="RT" & simdf$treat=="T"),]$period <- "2"

  simdf[which(simdf$sequence=="TR" & simdf$treat=="T"),]$period <- "1"
  simdf[which(simdf$sequence=="TR" & simdf$treat=="R"),]$period <- "2"

  simdf <- simdf[order(as.numeric(simdf[,1])),]
  simdf$cmax <- as.numeric(simdf$cmax)
  

  return(simdf)
}
