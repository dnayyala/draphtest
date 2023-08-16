dir.lkhd      <- function(x, mu, alpha) {
    #' @param x Data matrix of size $p \times n$ generated from [p]-dimensional Dirichlet distribution. 
    #' @param mu A scalar parameter which represents the dispersion parameter of the Dirichlet distribution. 
    #' @param alpha A $p - 1 \times 1$ vector which is the mean parameter of Dirichlet distribution. 
    #' @return a log-likelihood function of Dirichlet distribution. 
    
    ## Check if data and parameter have the same dimension
    if (nrow(x) != length(alpha)){
        stop("Data and parameter dimensions do not match")
    }
    
    theta <- exp(alpha)/sum(exp(alpha))
    n     <- ncol(x)   # sample size of date
    lkhd  <- n*lgamma(mu) -n*sum(lgamma(mu*theta)) + sum(t(mu*theta-1)%*%log(x))
    
    return(lkhd)
}       
