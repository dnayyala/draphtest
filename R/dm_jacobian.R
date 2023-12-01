dm.jacobian <- function(x, mu, alpha, param = "alpha"){
    #' Jacobian of Dirichlet distribution with respect to parameter alpha
    #' @param x A matrix size of (p x n) which is generated from [p]-dimensional Dirichlet distribution. 
    #' @param mu Scalar parameter of the Dirichlet distribution which represents the dispersion of the data. 
    #' @param alpha Vector of length p which represents the log-transformed mean vector of Dirichlet distribution. 
    #' @param param Character string identifying the parameter for which the Jacobian is being computed. The possible values are "mu" or "alpha".
    #' @return If ```param``` is set as ```alpha```, the function will return a $(p-1) \times 1$ vector of Jacobian, else if ```param``` is set as ```mu```, the function will return a scalar of Jacobian.
    
    n  <- ncol(x)   # sample size of date
    p  <- nrow(x)
    
    theta <- exp(alpha - max(alpha))/sum(exp(alpha - max(alpha)))
    X.plus  <- colSums(x)
    
    ## Check to confirm data and parameter vector are of the same dimension
    if (length(alpha) != p){
        stop("Data and parameter dimensions do not match.")
    }

    jacobian <- numeric(nrow(x))
    
    
    if (param == "alpha"){
        #' Jacobian in terms of parameter alpha
        Q       <- matrix(theta, nrow=p, ncol=p, byrow=TRUE)
        diag(Q) <- diag(Q)-1
        Q       <- Q[-1,]
        m       <- rowSums(digamma(x + mu*theta)) - n*digamma(mu*theta)
        
        jacobian <- mu*theta[-1] * (Q%*%m)
        
    } else if (param == "mu"){
        #' Jacobian in terms of parameter mu
        jacobian <- sum(digamma(mu)-digamma(X.plus+mu) + colSums(theta* digamma(x + mu*theta) - theta*digamma(mu*theta)))
    }
    
    return(jacobian);
}
