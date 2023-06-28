#' ver0: Reparameterized Dirichlet distribution

dir.lkhd      <- function(x, mu, theta) {
    #' @param x is size of (p x n) which is generated from [p]-dimensional Dirichlet distribution. 
    #' @param mu is a constant term. 
    #' @param alpha is the parameter vector of Dirichlet distribution. 
    
    theta <- exp(alpha)/sum(exp(alpha))
    n     <- ncol(x)   # sample size of date
    lkhd  <- n*lgamma(mu) -n*sum(lgamma(mu*theta)) + sum(t(mu*theta-1)%*%log(x))
    
    return(lkhd)
}

