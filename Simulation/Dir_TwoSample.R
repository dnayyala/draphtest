
# --------------------------------------------------------------------------------------------------
# Simulation for empirical type I error & power of two-sample test 
#
# Distribution: Dirichlet distribution
# Description: 
#   STEP 1: Generate empirical null distribution and compute cut-off values
#   STEP 2: Compute empirical type I error (diff.rate = 0) or power (diff.rate > 0)
#
# Arguments (via commandArgs)
#   args[1]: p             - dimension of original data before projection (p >> k)
#   args[2]: n             - sample size
#   args[3]: k             - projection dimension (K < n)
#   args[4]: mu.x          - dispersion parameter for group X 
#   args[5]: mu.y          - dispersion parameter for group Y
#   args[6]: diff.rate     - effect size (Type I error = 0, Power = 0.05/0.1/0.2)
# --------------------------------------------------------------------------------------------------

rm(list=ls())
#setwd("~path")
args <- commandArgs(trailingOnly = TRUE)

# Load R packages
library(mvtnorm)
library(foreach)
library(doParallel)
library(doRNG)
library(gtools)
library(ICSNP)
library(draphtest)

## Set values
## args will be a vector of length three containing p, n, k, mu, and diff.rate
p            <- as.integer(args[1])
n            <- as.integer(args[2])
k            <- as.integer(args[3])
mu.x         <- as.integer(args[4])
mu.y         <- as.integer(args[5])
diff.rate    <- as.numeric(args[6])
type1err     <- 0.05   # Significance level

# Generate Dirichlet random variables
set.seed(1234)
alpha      <- c(0, rgamma(p-1, 3, 5))
theta.null <- exp(alpha - max(alpha))/sum(exp(alpha - max(alpha)))
theta.x    <- theta.null

# Generate parameter vector of theta in sample Y 
if (diff.rate == 0){
  theta.y = theta.null
} else {
  e       <-  min(alpha[-1])  # to give a difference between null & alternative
  n.equal <- p*(1-diff.rate) 
  alpha.y <- alpha + e*(1:p >= n.equal)
  theta.y <- exp(alpha.y - max(alpha.y))/sum(exp(alpha.y - max(alpha.y)))
}

x   <- t(rdirichlet(n, mu.x*theta.x))
y   <- t(rdirichlet(n, mu.y*theta.y))


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
avg.p.value = foreach(i=1:n.boots, .packages=c('gtools', 'ICSNP', 'draphtest'), .combine='rbind') %dorng%  {
        x.boot  <- t(rdirichlet(n, mu.x*theta.null))
        y.boot  <- t(rdirichlet(n, mu.y*theta.null))
  
  for(m in 1:m.random){
        RP.orth <- ortho.randproj(nrow=k, ncol=p, method = "norm", seed = NULL) # Random projection method using orthogonal for RAPTT 
        RP.prop <- random.proj(k, p) # The proposed random projection method for Wald and LRT
        
        rx.dir <- RP.prop %*% x.boot
        ry.dir <- RP.prop %*% y.boot
        
        rx.raptt  <- RP.orth %*% x.boot
        ry.raptt  <- RP.orth %*% y.boot
        
        r.null.theta <- RP.prop %*% theta.null
        r.null.dir   <- log(r.null.theta) - log(r.null.theta[1]);
        r.null.raptt <- RP.orth %*% theta.null
        
        # Estimates parameters      
        parm.x.est      <- dir.nr(rx.dir)
        alpha.x.est.dir <- unlist(parm.x.est$alpha)
        r.mu.x          <- parm.x.est$mu
        
        parm.y.est      <- dir.nr(ry.dir)
        alpha.y.est.dir <- unlist(parm.y.est$alpha)
        r.mu.y          <- parm.y.est$mu
        
        # Estimate common theta parameter 
        fit.common    <- dir.nr(x = rx.dir, y = ry.dir)
        alpha.est.dir <- unlist(fit.common$alpha)
        
        
        p.value.wald[m] <- dir.wald(x=rx.dir, y=ry.dir, mu.x = r.mu.x, mu.y = r.mu.y, 
                                        alpha.x = alpha.x.est.dir, alpha.y=alpha.y.est.dir, type="two")$p.value
        p.value.lrt[m] <- dir.lrt(x=rx.dir, y=ry.dir, mu.x=r.mu.x, mu.y= r.mu.y, alpha=alpha.est.dir, 
                                      alpha.x = alpha.x.est.dir, alpha.y=alpha.y.est.dir)$p.value
        p.value.raptt[m] <- HotellingsT2(t(rx.raptt), t(ry.raptt), test = "chi")$p.value
    
  }
        c(mean(p.value.wald), mean(p.value.lrt), mean(p.value.raptt)) 
    
}

cut.off <- apply(avg.p.value, 2, quantile, probs = type1err)

# STEP 2
p.value.wald  <- numeric(m.random)
p.value.lrt   <- numeric(m.random)
p.value.raptt <- numeric(m.random)

mean.p = foreach(i=1:n.total, .packages=c('gtools', 'ICSNP', 'draphtest'), .combine='rbind') %dorng%  {
    
    x.boot  <- t(rdirichlet(n, mu.x*theta.x))
    y.boot  <- t(rdirichlet(n, mu.y*theta.y))
    
  for(m in 1:m.random){
      RP.orth <- ortho.randproj(nrow=k, ncol=p, method = "norm", seed = NULL) # Random projection method using orthogonal for RAPTT 
      RP.prop <- random.proj(k, p) # The proposed random projection method for Wald and LRT
      
      rx.dir <- RP.prop %*% x.boot
      ry.dir <- RP.prop %*% y.boot
      
      rx.raptt  <- RP.orth %*% x.boot
      ry.raptt  <- RP.orth %*% y.boot
      
      r.null.theta <- RP.prop %*% theta.null
      r.null.dir   <- log(r.null.theta) - log(r.null.theta[1]);
      r.null.raptt <- RP.orth %*% theta.null
      
      #' Parameter Estimation     
      parm.x.est  <- dir.nr(rx.dir)
      alpha.x.est.dir <- unlist(parm.x.est$alpha)
      r.mu.x <- parm.x.est$mu
      
      parm.y.est  <- dir.nr(ry.dir)
      alpha.y.est.dir <- unlist(parm.y.est$alpha)
      r.mu.y <- parm.y.est$mu
      
      r.total <- cbind(rx.dir, ry.dir); 
      alpha.est.dir  <- unlist(dir.nr(r.total)$alpha)
     
      p.value.wald[m] <- dir.wald(x=rx.dir, y=ry.dir, mu.x = r.mu.x, mu.y = r.mu.y, 
                                  alpha.x = alpha.x.est.dir, alpha.y=alpha.y.est.dir, type="two")$p.value
      p.value.lrt[m] <- dir.lrt(x=rx.dir, y=ry.dir, mu.x=r.mu.x, mu.y= r.mu.y, alpha=alpha.est.dir, 
                                alpha.x = alpha.x.est.dir, alpha.y=alpha.y.est.dir)$p.value
      p.value.raptt[m] <- HotellingsT2(t(rx.raptt), t(ry.raptt), test = "chi")$p.value
  }
  c(mean(p.value.wald), mean(p.value.lrt), mean(p.value.raptt)) 
}


end.time <- Sys.time()
running.time <-  end.time - start.time; running.time

## WRIET A TABLE
wald.emp.result   <- mean(mean.p[,1] < cut.off[1]); wald.emp.result
lrt.emp.result    <- mean(mean.p[,2] < cut.off[2]); lrt.emp.result
raptt.emp.result  <- mean(mean.p[,3] < cut.off[3]); raptt.emp.result

emp.result <- cbind(wald.emp.result, lrt.emp.result, raptt.emp.result)
colnames(emp.result) <- c("Wald","LRT", "RAPTT"); emp.result

## Save the above table as csv file
write.csv(emp.result, file=paste0("Dir_TwoSample_result_p_", args[1], 
                                  "_n_",                     args[2], 
                                  "_k_",                     args[3], 
                                  "_mu.x_",                  args[4], 
                                  "_mu.y_",                  args[5], 
                                  "_DiffRate_",              args[6], ".csv"))
