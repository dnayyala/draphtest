dir.lrt <- function(x, y = NULL, mu.x, mu.y = NULL, alpha.x, alpha.y = NULL, alpha , type = "one"){
    
    p <- length(alpha) # the number of parameters
    
    if (is.null(mu.y)){
        type = "one"
    } else {
        type = "two"
    }

    if (type == "one"){
        ## Check if data and parameter dimensions match
        if (length(alpha.x) && length(alpha) != p){
            stop("Data and parameter dimensions do not match.")
        }
        
    #' test statistic of LRT 
    test.stat <- -2 * (
                 dir.lkhd(x=x, mu=mu.x, alpha = alpha) # log-likelihood under the null 
               - dir.lkhd(x=x, mu=mu.x, alpha = alpha.x)) # log-likelihood under the alternative
 
    p.value     <- pchisq(test.stat, df=p-1, lower.tail=FALSE) 
    
    } else if (type == "two"){
        
        if (length(alpha.x) && length(alpha.y) != p){
            stop("Data and parameter dimensions do not match.")
        }
        
    test.stat <-  -2* (dir.lkhd(x, mu.x, alpha) + dir.lkhd(y, mu.y, alpha)
                - dir.lkhd(x, mu.x, alpha.x)- dir.lkhd(y, mu.y, alpha.y))
    
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