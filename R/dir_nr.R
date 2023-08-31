dir.nr <- function(x, mu = NULL, alpha = NULL, tol = 1e-4, maxiters = 1e2){
    #' Newton-Rapshon algorithm to estimate parameter alpha (or theta)
    #' @parm x A data matrix of size $p \times n$ generated from [p]-dimensional Dirichlet distribution. 
    #' @param mu A scalar parameter which represents the dispersion parameter of the Dirichlet distribution
    #' @param alpha Vector of length p which represents the log-transformed mean vector of Dirichlet distribution. 
    #' @param tol Tolerance cutoff for determining convergence of the Newton-Raphson algorithm
    #' @param maxiters Maximum number of iterations for the algorithm. Default is set at 100.
    #' @param param Character string identifying the parameter for which the Jacobian is being computed. The possible values are "mu" or "alpha".
    #' @return the function will return a $p \times 1$ vector containing the estimate of the param alpha and a the estimate of mu
    
  
    p <- nrow(x)
    # Initiate the values         
    alpha1 <- list();

        if (is.null(alpha)){
            alpha1[[1]]     <- log(rowMeans(x)) - log(rowMeans(x)[1])
        } else {
            alpha1[[1]] <- alpha
        }

        ## Check if data and parameter dimensions agree
        if (length(alpha1[[1]]) != p){
            stop("Data and parameter dimensions do not match")
        }
        
        mu1 <- c()        
        # Initiate the values 
        if (is.null(mu)){
            theta  <- exp(alpha1[[1]])/sum(exp(alpha1[[1]]))        
            mu.hat <- (theta*(1-theta))/apply(x, 1, var)  - 1
            mu1[1] <-  mean(mu.hat)
        } else {
            mu1[1] <- mu
        }
 
        Dir.lkhd <- c();
        Dir.lkhd[1] <- dir.lkhd(x, mu=mu1[1], alpha=alpha1[[1]]) # loglikelihood of data x
        
        
        step.size.vector  <- c();
        n.step   <- 1; 
        eps      <- 1e2;   # Difference in function between iterations
        
        while (eps > tol & n.step <= maxiters){
            n.step    <- n.step + 1;
            step.size <- 1
            
            Q  <- solve(dir.hessian(x, mu=as.numeric(mu1[n.step-1]), alpha=alpha1[[n.step - 1]], param="alpha"), dir.jacobian(x, mu=as.numeric(mu1[n.step-1]), alpha=alpha1[[n.step - 1]], param="alpha"))
            alpha1[[n.step]] <-  c(0, alpha1[[n.step -1]][2:p] - step.size*Q);
            
            W  <- dir.jacobian(x, mu=as.numeric(mu1[n.step - 1]), alpha=alpha1[[n.step]], param = "mu")/dir.hessian(x, mu=as.numeric(mu1[n.step - 1]), alpha=alpha1[[n.step]], param = "mu")
            
            ## Check if updated estimate is positive for different step sizes.
            step.size.check <- 10^(-(1:6))[(sapply(1:6, \(x) as.numeric(mu1[n.step -1]) - (10^(-x))*W) >= 0)]
            if (length(step.size.check) == 0){
                step.size == 0
            } else {
                step.size <- step.size.check[1]
            }
            mu1[n.step] <- mu1[n.step - 1] - step.size*W            
            step.size.vector[n.step - 1] <- step.size
            
            Dir.lkhd[n.step] <- dir.lkhd(x, mu=as.numeric(mu1[n.step]), alpha=alpha1[[n.step]])
            eps <-  (Dir.lkhd[n.step] - Dir.lkhd[n.step - 1]) / sqrt(abs(Dir.lkhd[n.step -1 ]))
        }  
        # If last update is in the opposite direction, use penultimate update value             
        if(eps > 0 ){
            mu.est    <- mu1[n.step]
            alpha.est <- alpha1[n.step]
        } else {
            mu.est    <- mu1[n.step-1]
            alpha.est <- alpha1[n.step-1]
        }     
        return(return(list(mu=mu.est, alpha=alpha.est)))
}
