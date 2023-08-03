dir.jacobian <- function(x, mu, alpha, param = "alpha"){
    #' Jacobian of Dirichlet distribution with respect to parameter alpha
    #' @param x A matrix size of (p x n) which is generated from [p]-dimensional Dirichlet distribution. 
    #' @param mu Scalar parameter of the Dirichlet distribution which represents the dispersion of the data. 
    #' @param alpha Vector of length p which represents the log-transformed mean vector of Dirichlet distribution. 
    #" @param param Character string identifying the parameter for which the Jacobian is being computed. The possible values are "mu" or "alpha".
    #" @return a $p$ vector Jacobian with respect to mu or alpha 
    
    p          <- nrow(x)   # the number of parameters (theta
    ## Check to confirm data and parameter vector are of the same dimension
    if (length(alpha) != p){
        stop("Data and parameter dimensions do not match.")
    }
    n     <- ncol(x)   # sample size of x
    theta <- exp(alpha)/sum(exp(alpha))
    jacobian <- numeric(nrow(x))
    
 
    if (param == "alpha"){
        for(i in 1:p){
            jacobian[i] <- mu*sum(sum(log(x[i,]) - digamma(mu*theta[i]))* ifelse(i == i, theta[i]*(1-theta[i]), -theta[i]*theta[-i]))
        }
        
    } else if (param == "mu"){
        jacobian <- n*digamma(mu) + sum(theta*rowSums(log(x)) - n*theta*digamma(mu*theta))
    }
    
    return(jacobian);
}

