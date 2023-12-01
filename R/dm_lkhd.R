dm.lkhd <- function(x, mu, alpha) {
    #' @param x Data matrix of size $p \times n$ generated from [p]-dimensional Dirichlet-multinomial  distribution. 
    #' @param mu A scalar parameter which represents the dispersion parameter of the distribution. 
    #' @param alpha A $p - 1 \times 1$ vector which is the mean parameter of the distribution. 
    #' @return a scalar which is the log-likelihood of the Dirichlet-multinomial distribution. 
    
    ## Check if data and parameter have the same dimension
    if (nrow(x) != length(alpha)){
        stop("Data and parameter dimensions do not match")
    }
    
    n       <- ncol(x)   # sample size of date
    theta   <- exp(alpha - max(alpha))/sum(exp(alpha - max(alpha)))
    X.plus  <- colSums(x)
    lkhd    <- sum(lgamma(mu) + lgamma(X.plus+1) - lgamma(X.plus+mu) + colSums(lgamma(x + mu*theta) - lgamma(mu * theta) - lgamma(x+1)))
    
    return(lkhd)
}
