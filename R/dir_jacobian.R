dir.jacobian <- function(x, mu, alpha, param = "mu"){
    #' Jacobian of Dirichlet distribution with repect to paramater alpha
    #' @param x A matrix size of (p x n) which is generated from [p]-dimensional Dirichlet distribution. 
    #' @param mu Scalar parameter of the Dirichlet distribution which represents the dispersion of the data. 
    #' @param alpha Vector of length p which represents the log-transformed mean vector of Dirichlet distribution. 
    #" @param param Character string identifying the parameter for which the Jacobian is being computed. The possible values are "mu" or "alpha".
    
    #" @return Bich Na: What is the value returned by the function?
    
    p     <- nrow(x)   # the number of parameters (theta)
    ## Check to confirm data and parameter vector are of the same dimension
    if (length(alpha) != p){
        stop("Data and parameter dimensions do not match.")
    }
    n     <- ncol(x)   # sample size of x
    theta <- exp(alpha)/sum(exp(alpha))

    if (param == "alpha"){
        sec.term <- numeric(p)
        for (i in 1:p){
            sec.term[i] <- mu*theta[i]*sum(n*theta[-i]*digamma(mu*theta[-i]) - rowSums(x[-i,]))         
        }        
        jacobian <- mu*theta*(1-theta) *(-n*digamma(mu*theta) + rowSums(x)) + sec.term
    } else if (param == "mu"){
        jacobian <- n*digamma(mu) + sum(theta*rowSums(log(x)) - n*theta*digamma(mu*theta))
    }
    
    return(jacobian);
}


