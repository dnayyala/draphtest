dir.wald <- function(x, y, mu.x, mu.y = NULL, alpha = NULL, alpha.x, alpha.y = NULL,  type = "one"){
    #' @title Wald test statistic for Dirichlet distribution
    #' 
    #' @description Wald-type test to parameter alpha of the Dirichlet distribution
    #' 
    #' @param x A data matrix of size \eqn{p \times n} generated from \eqn{p}-dimensional Dirichlet distribution. 
    #' @param y A data matrix of size \eqn{p \times n} generated from \eqn{p}-dimensional Dirichlet distribution. 
    #' @param mu.x A scalar parameter of x which represents the dispersion parameter of the Dirichlet distribution
    #' @param mu.y A scalar parameter of y which represents the dispersion parameter of the Dirichlet distribution
    #' @param alpha Vector of length p which represents the log-transformed mean vector under the null hypothesis.
    #' @param alpha.x Vector of length p which represents the log-transformed mean vector of x.  
    #' @param alpha.y Vector of length p which represents the log-transformed mean vector of y. 
    #' @param type If ```type``` is  set as ```one```, the function will return a one sample test result, else if ```type``` is  set as ```two``` then the function will return two samples test result.
    #' @import stats
    #' @export
    #' @return the function will return test statistic, p value and the result of the Wald-type test.
    
  p   <- nrow(x) 
  n   <- ncol(x) 
 
  
  if (is.null(mu.y)){
      type == "one"
  } else {
      type == "two"
  }
  
  
  if (is.null(alpha)){
      type == "two"
  } else {
      type == "one"
  }
  
  if (type == "one"){
  ## Check if data and parameter dimensions match
  if (length(alpha.x) && length(alpha) != p){
      stop("Data and parameter dimensions do not match.")
  }
  
      fisher   <- -dir.hessian(x, mu, alpha.x, param = "alpha")/n  # Fisher information matrix of data x
      est.cov  <- solve(fisher)
      diff     <-  as.matrix(alpha.x[2:p] - alpha[2:p]) 
    
      test.stat   <- t(diff) %*% solve(est.cov) %*% diff
      p.value     <- pchisq(test.stat, df=p-1, lower.tail=FALSE) 
      
    } else if (type == "two"){

    if (length(alpha.x) && length(alpha.y) != p){
        stop("Data and parameter dimensions do not match.")
    }
        
      p           <- nrow(x) 
      n           <- ncol(x) 
      m           <- ncol(y) 
      
      fisher.x      <- - dir.hessian(x, mu.x, alpha.x, param = "alpha")/n # Fisher information matrix of data x
      fisher.y      <- - dir.hessian(y, mu.y, alpha.y, param = "alpha")/m # Fisher information matrix of data x
      
      est.cov.x     <-  solve(fisher.x)
      est.cov.y     <-  solve(fisher.y)
    
      diff        <-  as.matrix(alpha.x[2:p] - alpha.y[2:p]) 
      
      test.stat   <-   t(diff) %*% solve(est.cov.x + est.cov.y) %*% diff
      p.value     <- pchisq(test.stat, df=p-1, lower.tail=FALSE) 
    }  
  result <- NULL;
  
  if (p.value < 0.05) {
    result <- "Reject H0"
  } else {
    result <- "Fail to Reject H0"
  }
  return(list(test.statistic=test.stat, p.value=p.value, result=result))
}
