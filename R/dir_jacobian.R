dir.jacobian <- function(x, mu, alpha, param = "alpha"){
    #' Jacobian of Dirichlet distribution with respect to parameter alpha
    #' @param x A matrix size of (p x n) which is generated from [p]-dimensional Dirichlet distribution. 
    #' @param mu Scalar parameter of the Dirichlet distribution which represents the dispersion of the data. 
    #' @param alpha Vector of length p which represents the log-transformed mean vector of Dirichlet distribution. 
    #" @param param Character string identifying the parameter for which the Jacobian is being computed. The possible values are "mu" or "alpha".
    #" @return a Jacobian, the first-order partial derivatives with respect to param alpha, size of $(p-1)$ vector
    
    p  <- nrow(x)   # the number of parameters
    
    ## Check to confirm data and parameter vector are of the same dimension
    if (length(alpha) != p){
        stop("Data and parameter dimensions do not match.")
    }
    n     <- ncol(x)   # sample size of x
    theta <- exp(alpha)/sum(exp(alpha))
    jacobian <- numeric(nrow(x))
    
    
    if (param == "alpha"){
        
        Q       <- matrix(theta, nrow=p, ncol=p, byrow=TRUE)
        diag(Q) <- diag(Q)-1
        Q       <- Q[-1,]
        W       <- n*digamma(mu*theta)- rowSums(log(x))
        
        jacobian <- mu*theta[-1] * (Q%*%W)
    
        } else if (param == "mu"){
        jacobian <- n*digamma(mu) + sum(theta*rowSums(log(x)) - n*theta*digamma(mu*theta))
    }
    
    return(jacobian);
}

