dir.jacobian <- function(x, mu, alpha, param = "alpha"){
  #' @title Jacobian for Dirichlet distribution 
  #' 
  #' @description This function computes the Jacobian of Dirichlet distribution with respect to the two parameters: alpha and mu
  #' 
  #' @param x A matrix size of \eqn{p \times n} which is generated from \eqn{p}-dimensional Dirichlet distribution. 
  #' @param mu Scalar parameter of the Dirichlet distribution which represents the dispersion of the data. 
  #' @param alpha Vector of length p which represents the log-transformed mean vector of Dirichlet distribution. 
  #' @param param Character string identifying the parameter for which the Jacobian is being computed. The possible values are "mu" or "alpha".
  #' @import stats
  #' @export
  #' @return A Jacobian vector containing the first-order partial derivatives with respect to param alpha, of length \eqn{p-1}. For "mu", the function returns a scalar Jacobian.
  
  p  <- nrow(x)   # the number of parameters
  
  ## Check to confirm data and parameter vector are of the same dimension
  if (length(alpha) != p){
    stop("Data and parameter dimensions do not match.")
  }
  
  n     <- ncol(x)   # sample size of x
  theta <- exp(alpha - max(alpha))/sum(exp(alpha - max(alpha)))
  jacobian <- numeric(nrow(x))
  
  Q       <- matrix(theta, nrow=p, ncol=p, byrow=TRUE)
  diag(Q) <- diag(Q)-1
  Q       <- Q[-1,]
  W       <- rowSums(log(x)) - n*digamma(mu*theta) 
  
  if (param == "alpha"){
    jacobian <- -mu*theta[-1] * (Q%*%W)
    
  } else if (param == "mu"){
    jacobian <- n*digamma(mu) + theta%*%W
  }
  
  return(jacobian);
}
