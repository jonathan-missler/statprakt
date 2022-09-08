#fuglsang
#install.packages("tseries")
library(tseries)

data <- read.csv("fuglsang.csv", sep = ";", dec = ",")
colnames(data) <- c("subject", "sequence", "period", "treat", "cmax")

testmod <- lm(log(cmax)~treat + subject + sequence + period, data = data)
testmod2 <- lme(log(cmax) ~ treat + sequence + period, random = ~1|subject, data = data)
summary(testmod)
summary(testmod2)


mean(data[1:24,]$cmax)
sd(data[1:24,]$cmax)

#fig1 replication
buster_cmax(data$cmax[data$treat == "R"])

res <- buster_conf(data)

kpss.test(res[[3]])


mod <-lmer(log(cmax) ~ treat + sequence + period + (1|subject), data = data[1:38,])
summary(mod)

buster_confplot(res[[1]], res[[3]])

buster_mseplot(res[[2]])

buster_residplot(res[[4]][data$treat == "T"])




buster(data)





################

#canada
data <- read.table("can_data.txt", header = TRUE)
data <- data[,c(1,2,3,4,7)]
colnames(data) <- c("subject", "period", "treat", "sequence", "cmax")
data <- data[order(data$subject),]


#fig1 replication
buster_cmax(data$cmax[data$treat == "T"])

res <- buster_conf(data)

testres <- kpss.test(res[[3]])

buster_confplot(res[[1]], res[[3]])

buster_mseplot(res[[2]])

buster_residplot(res[[4]][data$treat == "T"])


lmemod <- lme(log(cmax) ~ treat + sequence + period, random = ~1|subject, data = data)

coef(lmemod)
data.frame(intervals(lmemod, which = "fixed")$fixed)["treatT",c(1,3)]

buster(data)



############################################################
#SATOWIB
satodata <- read.csv("satowib.csv", sep = "")

satomat <- t(as.matrix(satodata))
satomat <- satomat[2:nrow(satomat),2:ncol(satomat)]

similarities <- SaToWIB(satomat)
similarities <- similarities[order(similarities[,3]),]
colnames(similarities) <- c("Profile1", "Profile2", "Score32", "Ratio")

#Grundsätzlich gleiche Ergebnisse wie Fuglsang, nur leicht veränderte Reihenfolge, teils leicht 
#andere Werte, aber die gleichen paare werden erkannt. 

###########################################################
#SATOWIB canada 

library(dplyr)

candataprof <- read.table("candata_prof.txt", sep = "")

candataprof[candataprof == "BLQ"] <- 0
idxstring <- paste(candataprof$SubjectID, as.character(candataprof$Period), sep = "P")
row.names(candataprof) <- idxstring

candat2 <- candataprof[,5:ncol(candataprof)]
candat2 <- mutate_all(candat2, function(x) as.numeric(x))
canmat <- as.matrix(candat2)

cansato <- SaToWIB(canmat)


model <- lm(data$cmax ~ data$treat)
plot(model)


###################################################################
#SIMULATED DATA 

sims <- sim_data(36, 6.645, 0.2, -0.311, 1234, manipulation = FALSE)

simmod <- lm(log(cmax)~   treat   +sequence + period, data = sims)
summary(simmod)
plot(simmod)

confint(simmod, c("treatT"), level = 0.9)

buster(sims)

simsres <- buster_conf(sims)

buster(data)

buster_test(sims)

buster_test(data)

bla <- buster_conf(sims)

kpss.test(bla[[3]])

point_est <- bla[[3]]
point_est <- point_est[(ceiling(length(point_est)/4)+1):length(point_est)]

plot(point_est)

install.packages("strucchange")
library(strucchange)

ocus.res <- efp(bla[[3]] ~ 1, type = "OLS-CUSUM")
plot(ocus.res)
plot(ocus.res, alpha = 0.01, alt.boundary = TRUE)
## calculate corresponding test statistic
test <- sctest(ocus.res)



#############################################################################################
#tests

#--------------POINT ESTIMATES------------------
tests <- vector()

for(i in 1:50){
  simdat <- sim_data(36, 6.65, 0.2, -0.301, i+567, manipulation = FALSE)
  
  tests <- c(tests, buster_test(simdat, alpha = 0.05))
  
}
sum(tests)


tests <- vector()

for(i in 1:50){
  simdat <- sim_data(36, 6.65, 0.2, -0.301, i+567, manipulation = TRUE)
  
  tests <- c(tests, buster_test(simdat, alpha = 0.05))
  
}
sum(tests)

#------------KPSS-------------------
tests <- vector()

for(i in 1:50){
  simdat <- sim_data(36, 6.645, 0.2, -0.311, i+567, manipulation = FALSE)
  
  tests <- c(tests, kpss_test(simdat, alpha = 0.05))
  
}
sum(tests)


tests <- vector()

for(i in 1:50){
  simdat <- sim_data(36, 6.645, 0.2, -0.311, i+567, manipulation = TRUE)
  
  tests <- c(tests, kpss_test(simdat, alpha = 0.05))
  
}
sum(tests)

#---------------adf stationary--------------
tests <- vector()

for(i in 1:50){
  simdat <- sim_data(36, 6.645, 0.2, -0.311, i+567, manipulation = FALSE)
  
  tests <- c(tests, adf_test(simdat, alpha = 0.05))
  
}
sum(tests)


tests <- vector()

for(i in 1:50){
  simdat <- sim_data(36, 6.645, 0.2, -0.311, i+567, manipulation = TRUE)
  
  tests <- c(tests, adf_test(simdat, alpha = 0.05))
  
}
sum(tests)

#------------KPSS-------------------
tests <- vector()

for(i in 1:50){
  simdat <- sim_data(36, 6.645, 0.2, -0.311, i+567, manipulation = FALSE)
  
  tests <- c(tests, kpss_test(simdat, alpha = 0.05))
  
}
sum(tests)


tests <- vector()

for(i in 1:50){
  simdat <- sim_data(36, 6.645, 0.2, -0.311, i+567, manipulation = TRUE)
  
  tests <- c(tests, kpss_test(simdat, alpha = 0.05))
  
}
sum(tests)


#---------CMAX DIFF-------------------
tests <- vector()

for(i in 1:50){
  simdat <- sim_data(36, 6.645, 0.2, -0.311, i+567, manipulation = FALSE)
  
  tests <- c(tests, cmaxdiff_test(simdat, alpha = 0.05))
  
}
sum(tests)


tests <- vector()

for(i in 1:50){
  simdat <- sim_data(36, 6.645, 0.2, -0.311, i+567, manipulation = TRUE)
  
  tests <- c(tests, cmaxdiff_test(simdat, alpha = 0.05))
  
}
sum(tests)


#------------------STRUCTURAL CHANGE----------------------
tests <- vector()

for(i in 1:50){
  simdat <- sim_data(36, 6.645, 0.2, -0.311, i+567, manipulation = FALSE)
  
  tests <- c(tests, strucc_test(simdat, alpha = 0.05))
  
}
sum(tests)


tests <- vector()

for(i in 1:50){
  simdat <- sim_data(36, 6.645, 0.2, -0.311, i+567, manipulation = TRUE)
  
  tests <- c(tests, strucc_test(simdat, alpha = 0.05))
  
}
sum(tests)




cmaxdiff <- log(data$cmax[data$treat == "T"]) - mean(log(data$cmax[data$treat == "T"]))
cmaxdiff <- log(sims$cmax[sims$treat == "T"]) - mean(log(sims$cmax[sims$treat == "T"]))

plot(cmaxdiff)
kpss.test(cmaxdiff)
# Vergleichbare ergebnisse wie beim point estimate. KPSS test erkennt manipulation oft nicht, schlägt
#aber auch selten bei nichtmanipulierten daten aus.

#structural change test erzeugt bei n = 36 gute ergebnisse für alpha = 0.2, hat sonst aber gleiche
#probleme wie KPSS test (bei weniger n)


#-----------------------RULE BASED TEST------------------------
tests <- vector()

for(i in 1:50){
  simdat <- sim_data(36, 6.645, 0.2, -0.311, i+567, manipulation = FALSE)
  
  tests <- c(tests, confinttest(simdat, lag = 3))
  
}
sum(tests)


tests <- vector()

for(i in 1:50){
  simdat <- sim_data(36, 6.645, 0.2, -0.311, i+567, manipulation = TRUE)
  
  tests <- c(tests, confinttest(simdat,lag = 3))
  
}
sum(tests)

#Bester Test bisher, lag parameter von 3 in den meisten Fällen ausreichend. bei n<30 für lag = 3
#viele falsch negative ergebnisse, aber wenig falsch positive. Für lag=2, viele falsch positive, 
#kaum falsch negative ergebnisse.

#---------CONFIDENCE INTERVAL WIDTH-----------

tests <- vector()

for(i in 1:50){
  simdat <- sim_data(36, 6.645, 0.2, -0.311, i+567, manipulation = FALSE)
  
  tests <- c(tests, confwidth_test(simdat, alpha  = 0.1))
  
}
sum(tests)


tests <- vector()

for(i in 1:50){
  simdat <- sim_data(36, 6.645, 0.2, -0.311, i+567, manipulation = TRUE)
  
  tests <- c(tests, confwidth_test(simdat, alpha = 0.1))
  
}
sum(tests)

#Funktioniert überhaupt nicht, schlägt immer aus.

#------------------STRUCTURAL BREAK CONFIDENCE WIDTH--------------------

tests <- vector()

for(i in 1:50){
  simdat <- sim_data(36, 6.645, 0.2, -0.311, i+567, manipulation = FALSE)
  
  tests <- c(tests, confwidthstrucc_test(simdat, alpha  = 0.05))
  
}
sum(tests)


tests <- vector()

for(i in 1:50){
  simdat <- sim_data(36, 6.645, 0.2, -0.311, i+567, manipulation = TRUE)
  
  tests <- c(tests, confwidthstrucc_test(simdat, alpha = 0.05))
  
}
sum(tests)
# Funktioniert auch überhaupt nicht, gleiches Problem wie bei KPSS test über Breite der Konfidenzintervalle


#-------------adf test explosive------------------
tests <- vector()

for(i in 1:50){
  simdat <- sim_data(36, 6.645, 0.2, -0.311, i+567, manipulation = FALSE)
  
  tests <- c(tests, adfexp_test(simdat, alpha = 0.05))
  
}
sum(tests)


tests <- vector()

for(i in 1:50){
  simdat <- sim_data(36, 6.645, 0.2, -0.311, i+567, manipulation = TRUE)
  
  tests <- c(tests, adfexp_test(simdat, alpha = 0.05))
  
}
sum(tests)

exdat <- data[1:48,]

sd(log(exdat$cmax))
mean()

sd(log(exdat$cmax[exdat$treat == "T"]))
sd(log(exdat$cmax[exdat$treat == "R"]))

mean(log(exdat$cmax[exdat$treat == "T"]))
mean(log(exdat$cmax[exdat$treat == "R"]))


hist(log(sims$cmax[sims$treat == "T"]))
hist(log(sims$cmax[sims$treat == "R"]))

hist(log(sims$cmax))

sd(log(sims$cmax[sims$treat == "T"]))
sd(log(sims$cmax[sims$treat == "R"]))

sd(sims$cmax[sims$treat == "T"])
sd(sims$cmax[sims$treat == "R"])


#------------effect adf explosive simulation-----------------
tests <- vector()

for(i in 1:480){
  if(i %% 2 == 0){
    size = 48
  } else{
    size = 36
  }
  
  simdat <- sim_data(size, 6.49, 0.2, -0.311, i+567, manipulation = FALSE)
  
  tests <- c(tests, adfexp_test(simdat, alpha = 0.05))
  
}
sum(tests)


tests <- vector()

for(i in 1:20){
  if(i %% 2 == 0){
    size = 48
  } else{
    size = 36
  }
  
  simdat <- sim_data(size, 6.49, 0.2, -0.311, i+123, manipulation = TRUE)
  
  tests <- c(tests, adfexp_test(simdat, alpha = 0.05))
  
}
sum(tests)


ppv <- function(prev, se = 1, sp = 0.96){
  up <- prev * se
  lwr <- up+(1-prev)*(1-sp)
  return(up/lwr)
}

ppv(1/3)
ppv(0.1)
ppv(0.01)
ppv(0.001)




