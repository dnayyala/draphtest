dm.hessian <- function(x, mu, alpha, param = "alpha"){
  #' Hessian matrix of Dirichlet distribution with respect to parameter alpha
  #' @param x The data matrix of size $p \times n$ which is generated from [p]-dimensional Dirichlet distribution. 
  #' @param mu A positive scalar parameter which represents the dispersion parameter of the Dirichlet distribution.  
  #' @param alpha A $p \times 1$ parameter vector which is the log-transformed mean of Dirichlet distribution. 
  #' @param param A character string identifying the parameter for which the Jacobian is being computed. The possible values are "mu" or "alpha".
  #' @return If ```param``` is set as ```alpha```, the function will return a $(p-1) \times (p-1)$ Hessian matrix, else if ```param``` is set as ```mu```, the function will return a scalar of Hessian
  
  p <- nrow(x)   # the number of parameters (theta)
  
  ## Check if data and parameter dimensions match
  if (length(alpha) != p){
    stop("Data and parameter dimensions do not match.")
  }
  n <- ncol(x)   # sample size of x
  
  ## Convert the log-mean parameter to mean
  theta <- exp(alpha - max(alpha))/sum(exp(alpha - max(alpha)))
  hessian <- matrix(NA, nrow = p-1, ncol = p-1)
  if (param == "alpha"){

    m       <- rowSums(digamma(x + mu*theta)) - n*digamma(mu*theta)
    b       <- rowSums(trigamma(x + mu*theta)) - n*trigamma(mu*theta)
    Q       <- matrix(theta, nrow=p, ncol=p, byrow=TRUE)
    diag(Q) <- diag(Q)-1
    Q       <- Q[-1,]  
    
    #' Off-diagonal elements of Hessian Matrix
    for(u in 2:p){
      for(v in 2:p){
        e_uv <- replace(numeric(p), c(u,v), 1) 
        hessian[u-1,v-1] <- mu*theta[u]*theta[v]*((2*theta - e_uv)%*%m + mu*(t(theta*(theta-e_uv)) %*% b))
      }   
    }
    
    
    #' Diagonal elements of Hessian Matrix
    diag(hessian) <- -mu*theta[-1]*((1-2*theta[-1])*(Q%*%m) - mu*theta[-1]*(Q^2%*%b))
    
    
  } else if (param == "mu"){
    hessian    <- n*trigamma(mu) - sum(trigamma(X.plus + mu)) + (theta^2 %*% b)
  }
  
  return(hessian);
}

  
  