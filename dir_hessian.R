dir.hessian<- function(x, mu, alpha){
  #' Hessian matrix of Dirichlet distribution with respect to parameter alpha
  #' @x is size of (p x n) which is generated from [p]-dimensional Dirichlet distribution. 
  #' @param mu is constant term. 
  #' @param alpha is the parameter vector of Dirichlet distribution. 
  
  p <- nrow(x)   # the number of parameters (theta)
  n <- ncol(x)   # sample size of x
  
  theta   <- exp(alpha)/sum(exp(alpha))
  Hessian <- matrix(NA, nrow = p, ncol = p)  ## matrix(NA, nrow = p, ncol = p)
  
  #' Off-diagonal elements of Hessian Matrix
  for(i in 1:p){
      for(j in 1:p ){
            Hessian[i,j] <- mu*theta[i]*theta[j]*(
                n*(mu*theta[i]*(1-theta[i])*trigamma(mu*theta[i]) + digamma(mu*theta[i])*(1+2*theta[i])
                  +mu*theta[j]*(1-theta[j])*trigamma(mu*theta[j])+ digamma(mu*theta[j])*(1+2*theta[j])
                  -sum(theta[c(-i,-j)]^2*mu*trigamma(mu*theta[c(-i,-j)]) - 2*digamma(mu*theta[c(-i,-j)])))
                +sum(log(x[i,])*(2*theta[i]-1) + log(x[j,])*(2*theta[j]-1) +2*sum(theta[c(-i,-j)]*log(x[c(-i,-j),]))))
  } }
  
  #' Diagonal elements of Hessian Matrix
  for(i in 1:p){
      Hessian[i,i] <- mu*theta[i]*
          (n*((1-theta[i])*(digamma(mu*theta[i])*(1+mu*theta[i]) - mu*theta[i]*trigamma(mu*theta[i]))
              - sum(theta[-i]*(mu*trigamma(mu*theta[-i])*theta[i]*theta[-i]-digamma(mu*theta[-i]))))
           + sum((1-theta[i])*(1-2*theta[i])*log(x[i,]) - colSums(theta[-i]*log(x[-i,]))))
  }

  return(Hessian);
}


dir.mu.hessian <- function(x, mu, alpha){
  #' Hessian matrix of Dirichlet distribution with respect to constant term mu
  #' @x is size of (p x n) which is generated from [p]-dimensional Dirichlet distribution. 
  #' @param mu is constant term
  #' @param alpha is the parameter vector of Dirichlet distribution 
  
  
  p <- nrow(x)   # the number of parameters (theta)
  n <- ncol(x)
  theta <- exp(alpha)/sum(exp(alpha))
  Hessian.mu   <- matrix(NA, nrow = p-1, ncol = p-1)  ## matrix(NA, nrow = p, ncol = p)
  
  Hessian.mu    <- n*(trigamma(mu) - sum(theta^2 * trigamma(mu*theta) ))
  
  return(Hessian.mu);
}
