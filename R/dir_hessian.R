dir.hessian<- function(x, mu, alpha, param = "alpha"){
  #' Hessian matrix of Dirichlet distribution with respect to parameter alpha
  #' @param x The data matrix of size $p \times n$ which is generated from [p]-dimensional Dirichlet distribution. 
  #' @param mu A positive scalar parameter which represents the dispersion parameter of the Dirichlet distribution.  
  #' @param alpha A $p \times 1$ parameter vector which is the log-transformed mean of Dirichlet distribution. 
  #' @param param A character string identifying the parameter for which the Jacobian is being computed. The possible values are "mu" or "alpha".

  #" @return Hessian Bich Na: What is the value returned by the function?
  
  p <- nrow(x)   # the number of parameters (theta)
  ## CHeck if data and parameter dimensions match
  if (length(alpha) != p){
    stop("Data and parameter dimensions do not match.")
  }
  n <- ncol(x)   # sample size of x

  ## Convert the log-mean parameter to mean
  theta   <- exp(alpha)/sum(exp(alpha))

  if (param == "alpha"){
    Hessian <- matrix(NA, nrow = p, ncol = p)
    
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
  } else if (param == "mu"){
    Hessian   <- matrix(NA, nrow = p-1, ncol = p-1)  ## matrix(NA, nrow = p, ncol = p) 
    Hessian    <- n*(trigamma(mu) - sum(theta^2 * trigamma(mu*theta) ))
  }

  return(Hessian);
}
