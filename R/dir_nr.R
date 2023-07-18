dir.nr <- function(x, mu, tol){
    #' Newton-Rapshon algorithm to estimate parameter alpha (or theta)
    #' @x is size of (p x n) which is generated from [p]-dimensional Dirichlet distribution. 
    #' @mu is constant term. 
    #' @theta is the parameter vector of Dirichlet distribution. The sum of parameters is 1. 
    #'
    
    alpha1            <- list();
    step.size.vector  <- c();
    Dir.lkhd          <- c();
    p <- nrow(x)
    
    # Initiate the values 
    alpha1[[1]]     <- log(rowMeans(x)) -log(rowMeans(x)[1])
    Dir.lkhd[1]     <- dir.lkhd(x, mu, alpha=alpha1[[1]]) # loglikelihood of data x

 
    n.step   <- 1; 
    eps      <- 1e2;   # Difference in function between iterations
    
    while (eps > tol & n.step <= 100){
        n.step <- n.step + 1;
        step.size <- 1
        
        Q  <- solve(dir.hessian(x, mu=mu, alpha=alpha1[[n.step - 1]]), dir.jacobian(x, mu=mu, alpha=alpha1[[n.step - 1]]))
        alpha1[[n.step]] <-  alpha1[[n.step -1]] - step.size*Q;
        
        
        Dir.lkhd[n.step] <- dir.lkhd(x, mu, alpha=alpha1[[n.step]])
        eps <-  (Dir.lkhd[n.step] - Dir.lkhd[n.step - 1]) / sqrt(abs(Dir.lkhd[n.step -1 ]))
        # print(eps)
    }
    
    #create  function of Dir X. MLE should be output. Using NR method - make function
    
    if(eps > 0 ){
        alpha.est <- alpha1[[n.step]]
    } else {
        alpha.est <- alpha1[[n.step - 1]]
    }
    
    theta.est <- exp(alpha.est)/sum(exp(alpha.est))
    return(list(alpha.est=alpha.est, theta.est=theta.est))
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






