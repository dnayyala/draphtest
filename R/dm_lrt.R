dm.lrt <- function(x, y = NULL, mu.x, mu.y = NULL, alpha.x, alpha.y = NULL, alpha.null , type = "one"){
    #' @title Likelihood ratio test for Dirichlet multinomial distribution
    #' 
    #' @description Likelihood-ratio test to parameter alpha of the Dirichlet-multinomial distribution
 
    #' @param x A data matrix of size \eqn{p \times n} generated from [p]-dimensional Dirichlet-Multinomial distribution. 
    #' @param y A data matrix of size \eqn{p \times n} generated from [p]-dimensional Dirichlet-Multinomial distribution. 
    #' @param mu.x A scalar parameter of x which represents the dispersion parameter.
    #' @param mu.y A scalar parameter of y which represents the dispersion parameter.
    #' @param alpha.null Vector of length p which represents the log-transformed mean vector under the null hypothesis.
    #' @param alpha.x Vector of length p which represents the log-transformed mean vector of x.  
    #' @param alpha.y Vector of length p which represents the log-transformed mean vector of y. 
    #' @param type If ```type``` is  set as ```one```, the function will return a one sample test result, else if ```type``` is  set as ```two``` then the function will return two samples test result.
    #' @import stats
    #' @export
    #' @return the function will return test statistic, p value and the result of the Likelihood-ratio test.

    p <- length(alpha.null) # the number of parameters
    n <- ncol(x) # the number of parameters
    
    if (is.null(mu.y)){
        type = "one"
    } else {
        type = "two"
    }
    
    if (type == "one"){
        ## Check if data and parameter dimensions match
        if (nrow(x) && length(alpha.null) != p){
            stop("Data and parameter dimensions do not match.")
        }
 

    #' test statistic of LRT 
    test.stat <- -2 * (
                 dm.lkhd(x=x, mu=mu.x, alpha = alpha.null)  # log-likelihood under the null
               - dm.lkhd(x=x, mu=mu.x, alpha = alpha.x)) # log-likelihood under the alternative

     p.value     <- pchisq(test.stat, df=p-1, lower.tail=FALSE) 

    
    } else if (type == "two"){
    
    if (length(alpha.x) && length(alpha.y) != p){
        stop("Data and parameter dimensions do not match.")
    }
 
     test.stat <- -2* (dm.lkhd(x, mu.x, alpha.null) + dm.lkhd(y, mu.y, alpha.null)
            - dm.lkhd(x, mu.x, alpha.x)- dm.lkhd(y, mu.y, alpha.y))

     p.value <- pchisq(test.stat, df=p-1, lower.tail=FALSE)
    
     }
    
     result <- NULL;

     if (p.value < 0.05) {
       result <- "Reject H0"
     } else {
       result <- "Fail to Reject H0"
     }
 return(list(test.statistic=test.stat, p.value=p.value, result=result))
}
