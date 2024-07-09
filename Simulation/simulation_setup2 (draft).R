# Simulation Setup 2: Simulation for empirical type 1 error of One sample test (Dirichlet-Multinomial distribution)

####################################################
##############         LIBRARY      ################
####################################################
rm(list=ls())
setwd("/Users/bchoi/NewFolder/draphtest/R")


#args <- commandArgs(trailingOnly = TRUE)
## args will be a vector of length three containing p n and k

# R packages
library(mvtnorm)
library(foreach)
library(doParallel)
library(doRNG)
library(gtools)
library(ICSNP)
library(draphtest)

start.time <- Sys.time()
## Set values
p        <- 100 #as.integer(args[1])
n        <- 50 ##as.integer(args[2])    # the number of samples from Dirichlet distribution
k        <- 15 #as.integer(args[3])

# Generate Dirichlet random variables
set.seed(1234)
mu     <- 100
alpha  <- c(0, rgamma(p-1, 3, 2))
theta.null <- exp(alpha - max(alpha))/sum(exp(alpha - max(alpha)))

type1err = 0.05


# Generate Dirichlet-multinomial distribution
lambda <- 100*p  #1e4 # empirical type 1 error (alpha) is related to the lambda value
x <- matrix(NA, p, n)

X.plus <- rpois(n, lambda) # X+ in the manuscript
for(i in 1:n){
  prob <- rdirichlet(1, mu*theta.null)
  x[,i]  <- rmultinom(1, size = X.plus[i], prob)
}


## STAGE 1
n.boots  <- 1e3
m.random <- 1e3
n.total  <- 1e3

n.cores <- detectCores()*0.75
registerDoParallel(cores=n.cores)

# STEP 1
p.value.wald  <- numeric(m.random)
p.value.lrt   <- numeric(m.random)
p.value.raptt <- numeric(m.random)


fnc.adj <- function(x, X.plus){
  #' The adjusted function is to conduct a Hotelling's T2 test (Mean vector testing)

  n <- length(X.plus)
  output <- matrix(NA, nrow(x), n)
  for(i in 1:n){
    output[,i] <- x[,i]/X.plus[i]
  }
  return(output)
}

avg.p.value = foreach(i=1:n.boots, .packages=c('gtools', 'ICSNP'), .combine='rbind') %dorng%  {
  x.boot <- matrix(NA, p, n)
  for(i in 1:n){
    prob       <- rdirichlet(1, mu*theta.null)
    x.boot[,i] <- rmultinom(1, size = X.plus[i], prob)
  }
  
  #' Adjustment bootsratp sample for RAPTT
  x.raptt   <- fnc.adj(x.boot, X.plus) # convert the projected data to mean vector for Hotelling's T2 test

  for(m in 1:m.random){
    RP.orth   <- ortho.randproj(nrow=k, ncol=p, method = "norm", seed = NULL) # Random projection method using orthogonal for RAPTT
    RP.prop   <- random.proj(p, k) # The proposed random projection method for Wald and LRT
    
    rx.dirmult <- RP.prop %*% x.boot
    rx.raptt <- (RP.orth %*% x.raptt)

    r.null.theta   <- c(RP.prop %*% theta.null)
    r.null.dirmult <- log(r.null.theta) - log(r.null.theta[1]);
    r.null.raptt   <- RP.orth %*% theta.null
    
    ## estimate parameter with the projected data
    parm.est  <- dm.nr(rx.dirmult)
    alpha.est <- unlist(parm.est$alpha)
    r.mu.boot <- parm.est$mu
    
    p.value.wald[m]  <- dm.wald(x=rx.dirmult, mu.x = r.mu.boot, alpha.null = r.null.dirmult, alpha.x=alpha.est)$p.value
    p.value.lrt[m]   <- dm.lrt(x=rx.dirmult, mu.x = r.mu.boot, alpha.null = r.null.dirmult, alpha.x=alpha.est)$p.value
    p.value.raptt[m] <- HotellingsT2(t(rx.raptt), mu = r.null.raptt, test = "chi")$p.value
  }
  c(mean(p.value.wald), mean(p.value.lrt), mean(p.value.raptt))
}


## avg.p.value is a matrix of size n.boots x 3

func.cutoff <- function(x){
  quantile(x, type1err)
}
cut.off <- apply(avg.p.value, 2, func.cutoff)

# STEP 2
#emp.alpha <- numeric(length(mu.Y))
p.value <- numeric(m.random)
mean.p <- numeric(n.total)

mean.p = foreach(i=1:n.total, .packages=c('gtools', 'ICSNP'), .combine='rbind') %dorng%  {
  
  x.boot <- matrix(NA, p, n)
  for(i in 1:n){
    prob <- rdirichlet(1, mu*theta.null)
    x.boot[,i]           <- rmultinom(1, size = X.plus[i], prob)
  }
  
  #' Adjustment bootsratp sample for RAPTT
  x.raptt   <- fnc.adj(x.boot, X.plus) # convert the projected data to mean vector for Hotelling's T2 test
  
  for(m in 1:m.random){
    RP.orth   <- ortho.randproj(nrow=k, ncol=p, method = "norm", seed = NULL) # Random projection method using orthogonal for RAPTT
    RP.prop   <- random.proj(p, k) # The proposed random projection method for Wald and LRT
    
    rx.dirmult <- RP.prop %*% x.boot
    rx.raptt   <- (RP.orth %*% x.raptt)

    r.null.theta   <- c(RP.prop %*% theta.null)
    r.null.dirmult <- log(r.null.theta) - log(r.null.theta[1]);
    r.null.raptt   <- RP.orth %*% theta.null
    
    ## estimate parameter with the projected data
    parm.est  <- dm.nr(rx.dirmult)
    alpha.est <- unlist(parm.est$alpha)
    r.mu.boot <- parm.est$mu
    
    p.value.wald[m]  <- dm.wald(x=rx.dirmult, mu.x = r.mu.boot, alpha.null = r.null.dirmult, alpha.x=alpha.est)$p.value
    p.value.lrt[m]   <- dm.lrt(x=rx.dirmult, mu.x = r.mu.boot, alpha.null = r.null.dirmult, alpha.x=alpha.est)$p.value
    p.value.raptt[m] <- HotellingsT2(t(rx.raptt), mu = r.null.raptt, test = "chi")$p.value
  }
  c(mean(p.value.wald), mean(p.value.lrt), mean(p.value.raptt)) 
}

raptt.emp.alpha  <- mean(mean.p[,3] < cut.off[3]); raptt.emp.alpha
lrt.emp.alpha <- mean(mean.p[,2] < cut.off[2]); lrt.emp.alpha
wald.emp.alpha <- mean(mean.p[,1] < cut.off[1]); wald.emp.alpha

## WRITE THE SAVE FILE NAME INCLUDING THE VALUES PROVIDED IN ARGS
## SOMETHING LIKE "Result_args[1]_args[2]_args[3].RData"
#WRIET A TABLE
emp.alpha <- cbind(raptt.emp.alpha, lrt.emp.alpha, wald.emp.alpha)
colnames(emp.alpha)<- c("RAPTT", "LRT", "WALD")
emp.alpha

end.time <- Sys.time()

running.time <- end.time - start.time; running.time

## WRITE THE SAVE FILE NAME INCLUDING THE VALUES PROVIDED IN ARGS
## SOMETHING LIKE "Result_args[1]_args[2]_args[3].RData"
#WRIET A TABLE
emp.alpha <- cbind(raptt.emp.alpha, lrt.emp.alpha, wald.emp.alpha)
colnames(emp.alpha)<- c("RAPTT", "LRT", "WALD")


## WRITE THE SAVE FILE NAME INCLUDING THE VALUES PROVIDED IN ARGS
## SOMETHING LIKE "Result_args[1]_args[2]_args[3].RData"
#save.image(file=paste0("/Users/bchoi/Google Drive/Projects/Random Project/Final version (042024)/results/DM_OneSample_alpha_p_", args[1],  "_n_",  args[2], "_k_",  args[3],  "_mu_", args[4], ".RData"))
#write.csv(emp.alpha, file=paste0("/Users/bchoi/Google Drive/Projects/Random Project/Final version (042024)/results/DM_OneSample_alpha_p_", p,  "_n_",  n, "_k_",  k, "_mu_", mu, ".csv"))
