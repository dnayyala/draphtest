dir.jacobian <- function(x, mu, alpha){
    #' Jacobian of Dirichlet distribution with repect to paramater alpha
    #' @x is size of (p x n) which is generated from [p]-dimensional Dirichlet distribution. 
    #' @param mu is constant term. 
    #' @param alpha is the parameter vector of Dirichlet distribution. 
    
    p     <- nrow(x)   # the number of parameters (theta)
    n     <- ncol(x)   # sample size of x
    theta <- exp(alpha)/sum(exp(alpha))
    
    sec.term <- numeric()
    for(i in 1:p){
        sec.term[i] <- mu*theta[i]*sum(n*theta[-i]*digamma(mu*theta[-i]) - rowSums(x[-i,]))         
    }
    
    jacobian <- c();
    jacobian <- mu*theta*(1-theta) *(-n*digamma(mu*theta) + rowSums(x)) + sec.term
    
    return(jacobian);
}



dir.mu.jacobian <- function(x, mu, alpha){
  #' Jacobian of Dirichlet distribution with repect to constant term mu
  #' @x is size of (p x n) which is generated from [p]-dimensional Dirichlet distribution.
  #' @param mu is constant term.
  #' @param alpha is the parameter vector of Dirichlet distribution.

  
  p     <- nrow(x)   # the number of parameters (theta)
  n     <- ncol(x)   # sample size of x
  theta <- exp(alpha)/sum(exp(alpha))

  jacobian.mu <- n*digamma(mu) + sum(theta*rowSums(log(x)) - n*theta*digamma(mu*theta))

  return(jacobian.mu);
}

