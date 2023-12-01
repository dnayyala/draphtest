dm.wald <- function(x, y, mu.x, mu.y = NULL, alpha.null = NULL, alpha.x, alpha.y = NULL,  type = "one"){
    #' Wald-type test to parameter alpha
    #' @parm x A data matrix of size $p \times n$ generated from [p]-dimensional Dirichlet-Multinomial distribution. 
    #' @parm y A data matrix of size $p \times n$ generated from [p]-dimensional Dirichlet-Multinomial distribution. 
    #' @parm mu.x A scalar parameter of x which represents the dispersion parameter 
    #' @parm mu.y A scalar parameter of y which represents the dispersion parameter 
    #' @param alpha Vector of length p which represents the log-transformed mean vector under the null hypothesis.
    #' @param alpha.x Vector of length p which represents the log-transformed mean vector of x.  
    #' @param alpha.y Vector of length p which represents the log-transformed mean vector of y. 
    #' @param type If ```type``` is  set as ```one```, the function will return a one sample test result, else if ```type``` is  set as ```two``` then the function will return two samples test result.
    #' @return the function will return test statistic, p value and the result of the Wald-type test.
    
   
    p <- nrow(x) 
    n <- ncol(x) 
    X.plus <- colSums(x)
    N      <- mean(X.plus)  # Total sum 
    
    if (is.null(mu.y)){
        type == "one"
    } else {
        type == "two"
    }
    
    if (is.null(alpha.null)){
        type == "two"
    } else {
        type == "one"
    }
    
    if (type == "one"){
        ## Check if data and parameter dimensions match
        if (length(alpha.x) && length(alpha.null) != p){
            stop("Data and parameter dimensions do not match.")
        }
        
    fisher  <- -dm.hessian(x, mu, alpha.x, param = "alpha") # Fisher information matrix of data x 
    diff    <- as.matrix(alpha.x[2:p] - alpha.null[2:p]) 

    test.stat   <-  t(diff) %*% (fisher) %*% diff
    p.value     <-  pchisq(test.stat, df=p-1, lower.tail=FALSE) 
  
    } else if (type == "two"){
    
    m <- ncol(y) 
    
    y.plus <- colSums(x)
    M      <- mean(y.plus)  # Total sum 
    
    
    fisher.x  <- - dm.hessian(x, mu.x, alpha.x, param = "alpha")# Fisher information matrix of data x
    fisher.y  <- - dm.hessian(y, mu.y, alpha.y, param = "alpha")# Fisher information matrix of data x
    
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
