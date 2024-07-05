dir.hessian <- function(x, mu, alpha, param = "alpha"){
  #' Hessian matrix of Dirichlet distribution with respect to parameter alpha
  #' 
  #' @param x The data matrix of size $(p x n)$ which is generated from [p]-dimensional Dirichlet distribution. 
  #' @param mu A positive scalar parameter which represents the dispersion parameter of the Dirichlet distribution.  
  #' @param alpha A $(p x 1)$ parameter vector which is the log-transformed mean of Dirichlet distribution. 
  #' @param param A character string identifying the parameter for which the Jacobian is being computed. The possible values are "mu" or "alpha".
  #' @return the Hessian matrix of Dirichlet distribution size of $(p-1) x (p-1)$
  
  p <- nrow(x)   # the number of parameters (theta)
  
  #print(p) 
  
  ## Check if data and parameter dimensions match
  if (length(alpha) != p){
    stop("Data and parameter dimensions do not match.")
  }
  n <- ncol(x)   # sample size of x
  
  ## Convert the log-mean parameter to mean
  theta <- exp(alpha - max(alpha))/sum(exp(alpha - max(alpha)))
  hessian <- matrix(NA, nrow = p-1, ncol = p-1)
  if (param == "alpha"){
    
    W       <- rowSums(log(x)) - n*digamma(mu*theta)
    Q       <- matrix(theta, nrow=p, ncol=p, byrow=TRUE)
    diag(Q) <- diag(Q)-1
    Q       <- Q[-1,]  

    #' Off-diagonal elements of Hessian Matrix
    for(u in 2:p){
      for(v in 2:p){
        e_uv <- replace(numeric(p), c(u,v), 1) 
        hessian[u-1,v-1] <- mu*theta[u]*theta[v]*((2*theta - e_uv)%*%W - n*mu*(t(theta*(theta-e_uv)) %*% trigamma(mu*theta)))
      }
    }
    
    #' Diagonal elements of Hessian Matrix
    diag(hessian) <- -mu*theta[-1]*((1-2*theta[-1])*(Q%*%W) + n*mu*theta[-1]*(Q^2%*%trigamma(mu*theta)))

  } else if (param == "mu"){
    hessian    <- as.numeric(n*trigamma(mu) - n*(trigamma(mu*theta)%*%theta^2))
  }
  
  return(hessian);
}

