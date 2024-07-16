# Simulation Setup 5: Simulation for empirical power of One sample test (Dirichlet distribution)

####################################################
##############         LIBRARY      ################
####################################################
rm(list=ls())
setwd("/Users/bchoi/NewFolder/draphtest/R")

#setwd("~path")

#.rs.restartR() # restart R session
# args <- commandArgs(trailingOnly = TRUE)
# args will be a vector of length three containing p, n and k

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
p            <- 100
n            <- 50 # the number of samples from Dirichlet distribution
k            <- 10 #as.integer(args[1])
type1.err    <- 0.05   # Significance level
diff.rate    <- 0.05

set.seed(1234)
# Generate a Dirichlet distribution parameters
mu         <- 100 #as.integer(args[2])

alpha      <- c(0, rgamma(p-1, 3, 5)) 
theta.null <- exp(alpha - max(alpha))/sum(exp(alpha - max(alpha)))

if (diff.rate == 0){
  theta.alt = theta.null
} else {
  e    <- min(theta.null)/2 - min(theta.null)*0.1 # to add a difference between null & alternative
  
  n.equal  <- p*(1-diff.rate) 
  n.diff.u   <- ceiling((p*diff.rate/2))
  n.diff.l   <- round(p*diff.rate/2)
  
  theta.alt1 <- theta.null[1:n.equal]
  theta.alt2 <- theta.null[(n.equal+1) : (n.equal+n.diff.l)] - e
  theta.alt3 <- theta.null[(n.equal+n.diff.l+1) : (n.equal+n.diff.l+n.diff.u)] + e
  
  theta.alt = c(theta.alt1, theta.alt2, theta.alt3)
  
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


avg.p.value = foreach(i=1:n.boots, .packages=c('gtools', 'ICSNP'), .combine='rbind') %dorng%  {
  x.boot  <- t(rdirichlet(n, mu*theta.null))
  
  for(m in 1:m.random){
    RP.orth   <- ortho.randproj(nrow=k, ncol=p, method = "norm", seed = NULL) # Random projection method using orthogonal for RAPTT
    RP.prop   <- random.proj(p, k) # The proposed random projection method for Wald and LRT
    
    rx.dir    <- RP.prop %*% x.boot
    rx.raptt  <- RP.orth %*% x.boot
    
    r.null.theta <- RP.prop %*% theta.null
    r.null.dir   <- log(r.null.theta) - log(r.null.theta[1]);
    r.null.raptt  <- RP.orth %*% theta.null
    
    ## estimate parameter with the projected data
    parm.est   <- dir.nr(rx.dir)
    alpha.est  <- unlist(parm.est$alpha)
    r.mu.boot  <- parm.est$mu
    
    p.value.wald[m]  <- dir.wald(x=rx.dir, mu.x = r.mu.boot, alpha = r.null.dir, alpha.x=alpha.est)$p.value
    p.value.lrt[m]   <- dir.lrt(x=rx.dir, mu.x = r.mu.boot, alpha = r.null.dir, alpha.x=alpha.est)$p.value
    p.value.raptt[m] <- HotellingsT2(t(rx.raptt), mu = r.null.raptt, test = "chi")$p.value
  }
  c(mean(p.value.wald), mean(p.value.lrt), mean(p.value.raptt))
}

func.cutoff <- function(x){
  quantile(x, type1.err)
}

cut.off <- apply(avg.p.value, 2, func.cutoff)

# STEP 2
p.value <- numeric(m.random)
mean.p  <- numeric(n.total)

mean.p = foreach(i=1:n.total, .packages=c('gtools', 'ICSNP'), .combine='rbind') %dorng%  {
  
  x.boot  <- t(rdirichlet(n, mu*theta.alt)) # the Bootstrap sample under the alternative hypothesis
  
  for(m in 1:m.random){
    
    RP.orth   <- ortho.randproj(nrow=k, ncol=p, method = "norm", seed = NULL) # Random projection method using orthogonal for RAPTT
    RP.prop   <- random.proj(p, k) # The proposed random projection method for Wald and LRT
    
    rx.dir    <- RP.prop %*% x.boot
    rx.raptt  <- RP.orth %*% x.boot
    
    
    r.null.theta <- RP.prop %*% theta.null
    r.null.dir   <- log(r.null.theta) - log(r.null.theta[1]);
    r.null.raptt <- RP.orth %*% theta.null
    
    
    ## Estimate parameters with projected data
    parm.est   <- dir.nr(rx.dir)
    alpha.est  <- unlist(parm.est$alpha)
    r.mu.boot  <- parm.est$mu
    
    p.value.wald[m]  <- dir.wald(x=rx.dir, mu.x = r.mu.boot, alpha = r.null.dir, alpha.x=alpha.est)$p.value
    p.value.lrt[m]   <- dir.lrt(x=rx.dir, mu.x= r.mu.boot, alpha=r.null.dir, alpha.x=alpha.est)$p.value
    p.value.raptt[m] <- HotellingsT2(t(rx.raptt), mu = r.null.raptt, test = "chi")$p.value
  }
  c(mean(p.value.wald), mean(p.value.lrt), mean(p.value.raptt)) 
}


end.time <- Sys.time()

running.time <-  end.time - start.time; running.time
wald.emp.power   <- mean(mean.p[,1] < cut.off[1]); wald.emp.power
lrt.emp.power    <- mean(mean.p[,2] < cut.off[2]); lrt.emp.power
raptt.emp.power  <- mean(mean.p[,3] < cut.off[3]); raptt.emp.power

emp.power <- cbind(wald.emp.power, lrt.emp.power, raptt.emp.power)
colnames(emp.power) <- c("wald","lrt", "raptt"); emp.power
