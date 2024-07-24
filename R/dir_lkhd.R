dir.lkhd      <- function(x, mu, alpha) {
    #' @title Dirichlet likelihood function
    #' 
    #' @description Function to compute the Dirichlet likelihood function given a vector of observations and parameters of the distribution
  
    #' @param x Data matrix of size \eqn{p \times n} generated from \eqn{p}-dimensional Dirichlet distribution. 
    #' @param mu A scalar parameter which represents the dispersion parameter of the Dirichlet distribution. 
    #' @param alpha A vector of length \eqn{p-1} which is the mean parameter of Dirichlet distribution. 
    #' @import stats
    #' @export
    #' @return a log-likelihood function of Dirichlet distribution. 
    
    ## Check if data and parameter have the same dimension
    if (nrow(x) != length(alpha)){
        stop("Data and parameter dimensions do not match")
    }
    
    theta <- exp(alpha - max(alpha))/sum(exp(alpha - max(alpha)))
    n     <- ncol(x)   # sample size of date
    lkhd  <- n*lgamma(mu) -n*sum(lgamma(mu*theta)) + sum(t(mu*theta-1)%*%log(x))
    
    return(lkhd)
}       
