dm.nr <- function(x, mu = NULL, alpha = NULL, tol = 1e-4, maxiters = 1e2){
    #' Newton-Rapshon algorithm to estimate parameter alpha (or theta)
    #' @parm x A data matrix of size (p x n) generated from [p]-dimensional Dirichlet-multinomial distribution. 
    #' @param mu A scalar parameter which represents the dispersion parameter of the distribution
    #' @param alpha Vector of length p which represents the log-transformed mean vector of the distribution. 
    #' @param tol Tolerance cutoff for determining convergence of the Newton-Raphson algorithm
    #' @param maxiters Maximum number of iterations for the algorithm. Default is set at 100.
    #' @param param Character string identifying the parameter for which the Jacobian is being computed. The possible values are "mu" or "alpha".
    #' @return the function will return a (p x 1) vector containing the estimate of the param alpha and a the estimate of mu
    p      <- nrow(x)
    X.plus <- colSums(x)
    N      <- mean(X.plus)  # Total sum 
    
    # Initiate the values         
    alpha1 <- list();
    theta1      <- rowMeans(x)/N
    
    if (is.null(alpha)){
        alpha1[[1]] <- log(theta1) - log(theta1[1])
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
        var.x <- apply(x, 1, var)
        mu.hat <- (theta1*N*(N-theta1) - var.x)/(var.x - theta1*N*(1-theta1))
        mu1[1] <-  mean(mu.hat)
    } else {
        mu1[1] <- mu
    }
    
    DM.lkhd <- c();
    DM.lkhd[1] <- dm.lkhd(x, mu=mu1[1], alpha=alpha1[[1]]) 
    
    n.step   <- 1; 
    eps      <- 1e2;   # Difference in function between iterations
    step.size.vector  <- c();
    
    while (eps > tol & n.step <= maxiters){
        n.step    <- n.step + 1;
       
        ## Fix mu and update alpha first 
        Q  <- solve(dm.hessian(x, mu=as.numeric(mu1[n.step-1]), alpha=alpha1[[n.step-1]], param = "alpha"), dm.jacobian(x, mu=as.numeric(mu1[n.step-1]), alpha=alpha1[[n.step-1]], param = "alpha"))
        ## Check if updated estimate is increasing the likelihood and not decreasing.        
        step.size <- 1
        step.size.check <- TRUE
        while (step.size.check && step.size >= 1e-6){
            alpha1[[n.step]] <-  c(0, alpha1[[n.step -1]][2:p] - step.size*Q);
            
            dm.lkhd.test <- dm.lkhd(x, mu=as.numeric(mu1[n.step-1]), alpha=alpha1[[n.step]])
            if(is.nan(dm.lkhd.test))
            {
                alpha1[[n.step]] <- alpha1[[n.step - 1]]
            }
            step.size.check <- ((dm.lkhd(x, mu=as.numeric(mu1[n.step-1]), alpha=alpha1[[n.step]]) - dm.lkhd(x, mu=as.numeric(mu1[n.step-1]), alpha=alpha1[[n.step-1]])) <= 0)
            step.size <- step.size*1e-1
        }
        if (step.size.check){
            alpha1[[n.step]] <- alpha1[[n.step - 1]]
        }
        

        ## Fix alpha and update mu
        W  <- dm.jacobian(x, mu=as.numeric(mu1[n.step - 1]), alpha=alpha1[[n.step]], param = "mu")/dm.hessian(x, mu=as.numeric(mu1[n.step - 1]), alpha=alpha1[[n.step]], param = "mu")
        ## Check if updated estimate is positive for different step sizes.        
        step.size <- 1
        step.size.check <- 10^(-(1:6))[(sapply(1:6, \(x) as.numeric(mu1[n.step -1]) - (10^(-x))*W) >= 0)]
        if (length(step.size.check) == 0){
            step.size == 0
        } else {
            step.size <- step.size.check[1]
        }
        mu1[n.step] <- mu1[n.step - 1] - step.size*W  
        
        DM.lkhd[n.step] <- dm.lkhd(x, mu=as.numeric(mu1[n.step]), alpha=alpha1[[n.step]])
        eps <-  (DM.lkhd[n.step] - DM.lkhd[n.step - 1]) / sqrt(abs(DM.lkhd[n.step -1 ]))
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
