# Simulation Setup 8: Simulation for empirical power of two sample test (Dirichlet-multinomial distribution)

####################################################
##############         LIBRARY      ################
####################################################
rm(list=ls())
setwd("/Users/bchoi/NewFolder/draphtest/R")

# R packages
library(mvtnorm)
library(foreach)
library(doParallel)
library(doRNG)
library(gtools)
library(ICSNP)
library(draphtest)

## Set values
p        <- 100 
k        <- 5      #as.integer(args[1])
n        <- 50     # the number of samples from Dirichlet distribution
typ1err  <- 0.05   # Significance level
diff.rate <- 0.05  # the difference rate between null % alternative parameters 


# Generate Dirichlet-multinomial random variables
set.seed(1234)
mu.x      <- p
mu.y      <- p  #as.integer(args[3])

# Generate the parameters
alpha        <- c(0, rgamma(p-1, 3, 5))
theta.null   <- exp(alpha - max(alpha))/sum(exp(alpha - max(alpha)))


# Generate parameter vector of theta in sample Y 
if (diff.rate == 0){
  theta.y = theta.null
} else {
  e    <- min(theta.null)/2 - min(theta.null)*0.1 # to give a difference between null & alternative
  
  n.equal  <- p*(1-diff.rate) 
  n.diff.u   <- ceiling((p*diff.rate/2))
  n.diff.l   <- round(p*diff.rate/2)
  
  theta.alt1 <- theta.null[1:n.equal]
  theta.alt2 <- theta.null[(n.equal+1) : (n.equal+n.diff.l)] - e
  theta.alt3 <- theta.null[(n.equal+n.diff.l+1) : (n.equal+n.diff.l+n.diff.u)] + e
  
  theta.y = c(theta.alt1, theta.alt2, theta.alt3)
}

theta.x <- theta.null


lambda <- 100*p  #1e4; empirical type 1 error (alpha) is related to the lambda value
x <- matrix(NA, p, n)
y <- matrix(NA, p, n)

X.plus <- rpois(n, lambda) # X+ in the manuscript
for(i in 1:n){
  prob <- rdirichlet(1, mu.x*theta.x)
  x[,i]  <- rmultinom(1, size = X.plus[i], prob)
}

Y.plus <- rpois(n, lambda) # X+ in the manuscript
for(i in 1:n){
  prob <- rdirichlet(1, mu.y*theta.y)
  y[,i]  <- rmultinom(1, size = Y.plus[i], prob)
}

fnc.adj <- function(x, X.plus){
  #' The adjusted function is to conduct a Hotelling's T2 test (Mean vector testing)
  
  n <- length(X.plus)
  output <- matrix(NA, nrow(x), n)
  for(i in 1:n){
    output[,i] <- x[,i]/X.plus[i]
  }
  return(output)
}

start.time <- Sys.time()

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

avg.p.value  <- numeric(m.random)
avg.p.value = foreach(i=1:n.boots, .packages=c('gtools', 'ICSNP'), .combine='rbind') %dorng%  {
  
  #' Generate the bootstrap samples under the null hypothesis 
  x.boot <- matrix(NA, p, n)
  for(i in 1:n){
    prob       <- rdirichlet(1, mu.x*theta.null)
    x.boot[,i] <- rmultinom(1, size = X.plus[i], prob)
  }

  y.boot  <- matrix(NA, p, n)
  for(i in 1:n){
    prob       <- rdirichlet(1, mu.y*theta.null)
    y.boot[,i] <- rmultinom(1, size = Y.plus[i], prob)
  }
  
  #' Adjustment bootsratp sample for RAPTT
  #' Convert the projected data to mean vector for Hotelling's T2 test
  x.raptt   <- fnc.adj(x.boot, X.plus)
  y.raptt   <- fnc.adj(y.boot, Y.plus) 
  
  
  for(m in 1:m.random){
    RP.orth <- ortho.randproj(nrow=k, ncol=p, method = "norm", seed = NULL) # Random projection method using orthogonal for RAPTT 
    RP.prop <- random.proj(p, k) # The proposed random projection method for Wald and LRT
    
    rx.dir <- RP.prop %*% x.boot
    ry.dir <- RP.prop %*% y.boot
    
    rx.raptt  <- RP.orth %*% x.raptt
    ry.raptt  <- RP.orth %*% y.raptt
    
    r.null.theta <- RP.prop %*% theta.null
    r.null.dir   <- log(r.null.theta) - log(r.null.theta[1]);
    r.null.raptt <- RP.orth %*% theta.null
    
    #' Estimates parameters      
    parm.x.est  <- dm.nr(rx.dir)
    alpha.x.est.dir <- unlist(parm.x.est$alpha)
    r.mu.x <- parm.x.est$mu
    
    parm.y.est  <- dm.nr(ry.dir)
    alpha.y.est.dir <- unlist(parm.y.est$alpha)
    r.mu.y <- parm.y.est$mu
    
    r.total <- cbind(rx.dir, ry.dir); 
    alpha.est.dir  <- unlist(dm.nr(r.total)$alpha)
    
    p.value.wald[m] <- dm.wald(x=rx.dir, y=ry.dir, mu.x = r.mu.x, mu.y = r.mu.y, 
                                alpha.x = alpha.x.est.dir, alpha.y=alpha.y.est.dir, type="two")$p.value
    p.value.lrt[m] <- dm.lrt(x=rx.dir, y=ry.dir, mu.x=r.mu.x, mu.y= r.mu.y, alpha=alpha.est.dir, 
                              alpha.x = alpha.x.est.dir, alpha.y=alpha.y.est.dir)$p.value
    p.value.raptt[m] <- HotellingsT2(t(rx.raptt), t(ry.raptt), test = "chi")$p.value
    
  }
  c(mean(p.value.wald), mean(p.value.lrt), mean(p.value.raptt)) 
  
}

func.cutoff <- function(x){
  quantile(x, typ1err)
}

cut.off <- apply(avg.p.value, 2, func.cutoff)


# STEP 2
#emp.alpha <- numeric(length(mu.Y))
p.value.wald  <- numeric(m.random)
p.value.lrt   <- numeric(m.random)
p.value.raptt <- numeric(m.random)

mean.p = foreach(i=1:n.total, .packages=c('gtools', 'ICSNP'), .combine='rbind') %dorng%  {
  #' Generate the bootstrap samples under the alternative hypothesis
  #' 
  x.boot <- matrix(NA, p, n)
  for(i in 1:n){
    prob       <- rdirichlet(1, mu.x*theta.x)
    x.boot[,i] <- rmultinom(1, size = X.plus[i], prob)
  }
  
  y.boot  <- matrix(NA, p, n)
  for(i in 1:n){
    prob       <- rdirichlet(1, mu.y*theta.y)
    y.boot[,i] <- rmultinom(1, size = Y.plus[i], prob)
  }
  
  #' Adjustment bootsratp sample for RAPTT
  #' Convert the projected data to mean vector for Hotelling's T2 test
  x.raptt   <- fnc.adj(x.boot, X.plus)
  y.raptt   <- fnc.adj(y.boot, Y.plus) 
  
  for(m in 1:m.random){
    RP.orth <- ortho.randproj(nrow=k, ncol=p, method = "norm", seed = NULL) # Random projection method using orthogonal for RAPTT 
    RP.prop <- random.proj(p, k) # The proposed random projection method for Wald and LRT
    
    rx.dir <- RP.prop %*% x.boot
    ry.dir <- RP.prop %*% y.boot
    
    rx.raptt  <- RP.orth %*% x.raptt
    ry.raptt  <- RP.orth %*% y.raptt
    
    r.null.theta <- RP.prop %*% theta.null
    r.null.dir   <- log(r.null.theta) - log(r.null.theta[1]);
    r.null.raptt <- RP.orth %*% theta.null
    
    #' Parameter Estimation     
    parm.x.est  <- dm.nr(rx.dir)
    alpha.x.est.dir <- unlist(parm.x.est$alpha)
    r.mu.x <- parm.x.est$mu
    
    parm.y.est  <- dm.nr(ry.dir)
    alpha.y.est.dir <- unlist(parm.y.est$alpha)
    r.mu.y <- parm.y.est$mu
    
    r.total <- cbind(rx.dir, ry.dir); 
    alpha.est.dir  <- unlist(dm.nr(r.total)$alpha)
    
    p.value.wald[m] <- dm.wald(x=rx.dir, y=ry.dir, mu.x = r.mu.x, mu.y = r.mu.y, 
                                alpha.x = alpha.x.est.dir, alpha.y=alpha.y.est.dir, type="two")$p.value
    p.value.lrt[m] <- dm.lrt(x=rx.dir, y=ry.dir, mu.x=r.mu.x, mu.y= r.mu.y, alpha=alpha.est.dir, 
                              alpha.x = alpha.x.est.dir, alpha.y=alpha.y.est.dir)$p.value
    p.value.raptt[m] <- HotellingsT2(t(rx.raptt), t(ry.raptt), test = "chi")$p.value
  }
  c(mean(p.value.wald), mean(p.value.lrt), mean(p.value.raptt)) 
}


end.time <- Sys.time()
running.time <- end.time - start.time; running.time

wald.emp.power  <- mean(mean.p[,1] < cut.off[1]); wald.emp.power 
lrt.emp.power   <- mean(mean.p[,2] < cut.off[2]); lrt.emp.power
raptt.emp.power <- mean(mean.p[,3] < cut.off[3]); raptt.emp.power


emp.power <- cbind(wald.emp.power, lrt.emp.power, raptt.emp.power)
colnames(emp.power) <- c("WALD","LRT", "RAPTT"); emp.power
