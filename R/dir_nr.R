dir.nr <- function(x, mu = NULL, alpha = NULL, tol = 1e-4, , maxiters = 1e2, param = "alpha"){
    #' Newton-Rapshon algorithm to estimate parameter alpha (or theta)
    #' @parm x A data matrix of size $p \times n$ generated from [p]-dimensional Dirichlet distribution. 
    #' @param mu A scalar parameter which represents the dispersion parameter of the Dirichlet distribution
    #' @param alpha Vector of length p which represents the log-transformed mean vector of Dirichlet distribution. 
    #' @param tol Tolerance cutoff for determining convergence of the Newton-Raphson algorithm
    #' @param maxiters Maximum number of iterations for the algorithm. Default is set at 100.
    #' @param param Character string identifying the parameter for which the Jacobian is being computed. The possible values are "mu" or "alpha".

    #' @return If ```param``` is set as ```alpha```, the function will return a $p \times 1$ vector containing the estimate of the 
    
    step.size.vector  <- c();
    Dir.lkhd          <- c();
    p <- nrow(x)
    ## Check if data and parameter dimensions agree
    if (length(alpha) != p){
        stop("Data and parameter dimensions do not match")
    }
 
    n.step   <- 1; 
    eps      <- 1e2;   # Difference in function between iterations

    if (param == "alpha"){
        # Initiate the values         
        alpha1            <- list();
        alpha1[[1]]     <- log(rowMeans(x)) - log(rowMeans(x)[1])
        Dir.lkhd[1]     <- dir.lkhd(x, mu, alpha=alpha1[[1]]) # loglikelihood of data x
        
        while (eps > tol & n.step <= maxiters){
            n.step <- n.step + 1;
            step.size <- 1
            
            Q  <- solve(dir.hessian(x, mu=mu, alpha=alpha1[[n.step - 1]], param), dir.jacobian(x, mu=mu, alpha=alpha1[[n.step - 1]], param))
            alpha1[[n.step]] <-  alpha1[[n.step -1]] - step.size*Q;
            
            
            Dir.lkhd[n.step] <- dir.lkhd(x, mu, alpha=alpha1[[n.step]])
            eps <-  (Dir.lkhd[n.step] - Dir.lkhd[n.step - 1]) / sqrt(abs(Dir.lkhd[n.step -1 ]))
            # print(eps)
        }
        
        # If last update is in the opposite direction, use penultimate update value     
        if (eps > 0 ){
            alpha.est <- alpha1[[n.step]]
        } else {
            alpha.est <- alpha1[[n.step - 1]]
        }

        ## Return the final NR estiamte for alpha
        return(alpha.est)
    } else if (param == "mu"){
        mu1 <- list()        
        # Initiate the values         
        theta <- exp(alpha)/sum(exp(alpha))        
        mu.hat <- (theta*(1 - theta))/apply(x, 1, var)  - 1
        mu1[1]          <-  mean(mu.hat) 
        Dir.lkhd[1]     <- dir.lkhd(x, mu1[1], alpha) # loglikelihood of data x        
        while (eps > tol & n.step <= maxiters){   
            n.step <- n.step + 1;
            Q  <- dir.jacobian(x, mu=mu1[n.step - 1], alpha, param = "mu")/dir.mu.hessian(x, mu=mu1[n.step - 1], alpha, param = "mu")

            ## Check if updated estimate is positive for different step sizes.
            step.size.check <- 10^(-(1:6))[(sapply(1:6, \(x) mu1[n.step -1] - (10^(-x))*Q) >= 0)]
            if (length(step.size.check) == 0){
                step.size == 0
            } else {
                step.size <- step.size.check[1]
            }
            mu1[n.step] <- mu[n.step - 1] - step.size*Q            
            step.size.vector[n.step - 1] <- step.size
            
            Dir.lkhd[n.step] <- dir.lkhd(x, mu = mu1[n.step], alpha = alpha)
            eps <- (Dir.lkhd[n.step] - Dir.lkhd[n.step - 1]) / sqrt(abs(Dir.lkhd[n.step -1 ]))
            # print(eps)
        }
    
        # If last update is in the opposite direction, use penultimate update value             
        if(eps > 0 ){
            mu.est <- mu1[n.step]
        } else {
            mu.est <- mu1[n.step-1]
        }
        return(mu.est)
    }
}




dir.mu <- function(x, alpha, tol){
    #' Newton-Rapshon algorithm to estimate constant term mu.
    #' @x is size of (p x n) which is generated from [p]-dimensional Dirichlet distribution. 
    #' @mu is constant term. 
    #' @theta is the parameter vector of Dirichlet distribution. The sum of parameters is 1. 
    
    
    mu1               <- c();
    step.size.vector  <- c();
    Dir.lkhd          <- c();
    p <- nrow(x)
    theta <- exp(alpha)/sum(exp(alpha))
    # Initiate the values 
    
    variance <-apply(x, 1, var)
    
    mu.hat <- c()
    for(i in 1:length(theta)){
        mu.hat[i] <- theta[i]*(1-theta[i])/variance[i] -1 # MME of mu
    }
    
    
    mu1[1]          <-  mean(mu.hat) 
    Dir.lkhd[1]     <- dir.lkhd(x, mu1[1], theta=theta) # loglikelihood of data x
    
    
    n.step   <- 1; 
    eps      <- 1e2;   # Difference in function between iterations
    
    
    while (eps > tol & n.step <= 100){
        step.size <- 1
        n.step <- n.step + 1;
        Q  <- dir.mu.jacobian(x, mu=mu1[n.step - 1], theta)/dir.mu.hessian(x, mu=mu1[n.step - 1], theta=theta)
        
        out.of.range <- TRUE;
        
        while (out.of.range){
            mu1[n.step] <-  mu1[n.step -1] - step.size*Q;
            out.of.range <- ((mu1[n.step] <= 0 ) != 0)
            step.size <- step.size*0.1;
        }
        
        step.size.vector[n.step - 1] <- step.size
        
        Dir.lkhd[n.step] <- dir.lkhd(x, mu1[n.step], theta=theta)
        eps <- (Dir.lkhd[n.step] - Dir.lkhd[n.step - 1]) / sqrt(abs(Dir.lkhd[n.step -1 ]))
        # print(eps)
    }
    
    #create  function of Dir X. MLE should be output. Using NR method - make function
    
    if(eps > 0 ){
        mu.est <- mu1[n.step]
    } else {
        mu.est <- mu1[n.step-1]
    }
    
    return(mu.est)
}






