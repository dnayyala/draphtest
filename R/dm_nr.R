dm.nr <- function(x, y = NULL, mu = NULL, mu.y = NULL, alpha = NULL, 
                  tol = 1e-4, maxiters = 1e2, delta = 1e-4){
  #' @title Newton-Raphson for Dirichlet-multinomial distribution
  #' 
  #' @description Newton-Rapshon algorithm to estimate parameters of the
  #' Dirichlet-multinomial distribution. Supports two modes:
  #' \itemize{
  #'  \item \strong{One-sample} (y = NULL): estimates mu and alpha jointly from x.
  #'  \item \strong{Two-sample} (y supplied): estimates mu.x and mu.y separately
  #'    from each group, then fixes them while estimating a common alpha
  #'    from the combined likelihood. This avoids bias when mu_x != mu.y.
  #'  }
  #' 
  #' @param x A data matrix of size \eqn{p \times n} generated from [p]-dimensional 
  #' Dirichlet-multinomial distribution for group X. 
  #' @param y A data matrix of size \eqn{p \times n} generated from [p]-dimensional 
  #' Dirichlet-multinomial distribution for gropu Y. 
  #' @param mu A scalar parameter which represents the dispersion parameter of the distribution
  #' @param mu.y Dispersion parameter for y. If NULL, estimated from data.
  #' @param alpha Vector of length p which represents the log-transformed mean vector of the distribution. 
  #' @param tol Tolerance cutoff for determining convergence of the Newton-Raphson algorithm
  #' @param maxiters Maximum number of iterations for the algorithm. Default is set at 100.
  #' @param delta Threshold below which theta components are held fixed.
  #' @import stats
  #' @export
  #' @return the function will return a vector of length \eqn{p} containing the estimate of the param alpha and a the estimate of mu
  
  
  
  #--------------------------------------------------------------------------
  # TWO-SAMPLE CASE
  # mu.x and mu.y estimated separately, then fixed to estimate common alpha
  #--------------------------------------------------------------------------
  
  
  if (!is.null(y)){
    # Step 1. Estimate mu.x and mu.y separately 
    fit.x    <- dm.nr(x, mu = mu, alpha = alpha, tol = tol, 
                      maxiters = maxiters, delta = delta)
    mu.x.est <- fit.x$mu
    
    fit.y    <- dm.nr(y, mu = mu.y, alpha = alpha, tol = tol, 
                      maxiters = maxiters, delta = delta)
    mu.y.est <- fit.y$mu
    
    # Step 2. Estimate common theta with the estimated mu.x and mu.y
    p      <- nrow(x)
    alpha1 <- list()
    
    if (is.null(alpha)){
      X.plus   <- mean(colSums(x))
      Y.plus   <- mean(colSums(y))
      theta_xy <- (rowMeans(x)/X.plus + rowMeans(y)/Y.plus)/2
      theta_xy <- pmax(theta_xy, 1e-10)
      alpha1[[1]]   <- log(theta_xy)- log(theta_xy[1])
    } else {
      alpha1[[1]] <- alpha
    }
    
    ## Check if data and parameter dimensions agree
    if (length(alpha1[[1]]) != p){
      stop("Data and parameter dimensions do not match")
    }
    
    
    # Combined likelihood helpers (mu fixed)
    lkhd.comb <- function(a)
      dm.lkhd(x, mu.x.est, a) + dm.lkhd(y, mu.y.est, a) 
    
    jcb.comb <- function(a)
      dm.jacobian(x, mu.x.est, a, param = "alpha") + 
      dm.jacobian(y, mu.y.est, a, param = "alpha") 
    
    hess.comb <- function(a)
      dm.hessian(x, mu.x.est, a, param = "alpha") + 
      dm.hessian(y, mu.y.est, a, param = "alpha") 
    
    lkhd.val <- lkhd.comb(alpha1[[1]])
    n.step   <- 1
    eps      <- 1e2
    
    while (eps > tol & n.step <= maxiters){
      n.step    <- n.step + 1;
      
      # Current theta vector
      theta_cur <- exp(alpha1[[n.step - 1]] - max(alpha1[[n.step - 1]]))
      theta_cur <- theta_cur/sum(theta_cur)
      fix_index <- which(theta_cur < delta)
      
      
      # Compute Hessian and Jacobian for alpha
      H <- hess.comb(alpha1[[n.step - 1]])
      J <- jcb.comb(alpha1[[n.step - 1]])
      
      # If the Q is singular, then return Q as 0, then do not update alpha 
      Q <- tryCatch(
        solve(H, J),
        error=function(e){
          warning("Hessian matrix is singular. Skip alpha update at step:", n.step)
          rep(0, length(J))
        }
      )
      
      # Line Search 
      # Check if updated estimate is increasing the likelihood     
      step.size <- 1
      step.size.check <- TRUE
      while (step.size.check && step.size >= 1e-6){
        alpha_cand      <- alpha1[[n.step -1]]
        alpha_cand[1]   <- 0
        alpha_cand[2:p] <- alpha1[[n.step -1]][2:p] - step.size*Q
        improvement     <- lkhd.comb(alpha_cand) - lkhd.comb(alpha1[[n.step -1]]) 
        step.size.check <- is.na(improvement) || improvement <= 0
        step.size       <- step.size*1e-1
      }
      
      alpha1[[n.step]] <- if (step.size.check) alpha1[[n.step-1]] else alpha_cand
      
      
      # Fix alpha if alpha has near-zero values
      if (length(fix_index)>0){
        alpha1[[n.step]][fix_index] <- alpha1[[n.step - 1]][fix_index]
      }
      
      lkhd.new <- lkhd.comb(alpha1[[n.step]])
      eps      <- (lkhd.new - lkhd.val)/sqrt(abs(lkhd.val))  
      lkhd.val <- lkhd.new  
    }
    
    alpha.est <- if (eps > 0) alpha1[[n.step]] else alpha1[[n.step-1]]
    
    return(list(mu=list(mu.x=mu.x.est, mu.y=mu.y.est), alpha=alpha.est))
  }
  
  #--------------------------------------------------------------------------
  # One-SAMPLE CASE
  #--------------------------------------------------------------------------
  
  p      <- nrow(x)
  X.plus <- colSums(x)
  N      <- mean(X.plus)  # Total sum 
  
  # Initiate the values         
  alpha1 <- list()
  theta1 <- rowMeans(x)/N
  theta1 <- pmax(theta1, 1e-10)
  
  if (is.null(alpha)){
    alpha1[[1]] <- log(theta1) - log(theta1[1])
  } else {
    alpha1[[1]] <- alpha
  }
  
  ## Check if data and parameter dimensions agree
  if (length(alpha1[[1]]) != p){
    stop("Data and parameter dimensions do not match")
  }
  
  # Initiate a dispersion parameter mu
  mu1 <- c()        
  if (is.null(mu)){
    var.x  <- apply(x, 1, var)
    var.x  <- pmax(var.x, 1e-10)
    
    mu.hat <- (theta1*N*(N-theta1) - var.x)/(var.x - theta1*N*(1-theta1))
    mu.hat <- mu.hat[is.finite(mu.hat)]
    mu1[1] <- max(mean(mu.hat), 1e-3) # prevent mu.hat goes negative value
  } else {
    mu1[1] <- mu
  }
  
  DM.lkhd    <- c();
  DM.lkhd[1] <- dm.lkhd(x, mu=mu1[1], alpha=alpha1[[1]]) 
  
  n.step   <- 1; 
  eps      <- 1e2;   # Difference in function between iterations
  step.size.vector  <- c();
  
  while (eps > tol & n.step <= maxiters){
    n.step    <- n.step + 1;
    
    
    # Sparse componenent check
    theta_cur <- exp(alpha1[[n.step - 1]] - max(alpha1[[n.step - 1]]))
    theta_cur <- theta_cur/sum(theta_cur)
    fix_index <- which(theta_cur < delta)
    
    # Compute Hessian and Jacobian for alpha
    H <- dm.hessian(x, mu=as.numeric(mu1[n.step-1]), 
                    alpha=alpha1[[n.step-1]], param = "alpha")
    J <- dm.jacobian(x, mu=as.numeric(mu1[n.step-1]), 
                     alpha=alpha1[[n.step-1]], param = "alpha")
    
    # If the Q is singular, then return Q as 0, then do not update alpha 
    Q <- tryCatch(
      solve(H, J),
      error=function(e){
        warning("Hessian matrix is singular. Skip alpha update at step:", n.step)
        rep(0, length(J))
      }
    )
    
    ## Check if updated estimate is increasing the likelihood and not decreasing.        
    step.size <- 1
    step.size.check <- TRUE
    while (step.size.check && step.size >= 1e-6){
      alpha1[[n.step]] <-  c(0, alpha1[[n.step -1]][2:p] - step.size*Q);
      
      dm.lkhd.test <- dm.lkhd(x, mu=as.numeric(mu1[n.step-1]), 
                              alpha=alpha1[[n.step]])
      if(is.nan(dm.lkhd.test)){
        alpha1[[n.step]] <- alpha1[[n.step - 1]]
        break
      }
      
      improvement     <- dm.lkhd.test - dm.lkhd(x, mu=as.numeric(mu1[n.step-1]), 
                                                alpha=alpha1[[n.step-1]])
      step.size.check <-  is.na(improvement) || improvement <= 0
      step.size       <- step.size * 1e-1
    }
    if (step.size.check){
      alpha1[[n.step]] <- alpha1[[n.step - 1]]
    }
    
    # Fix sparse components
    if (length(fix_index) > 0)
      alpha1[[n.step]][fix_index] <- alpha1[[n.step - 1]][fix_index]
    
    
    ## Update mu
    J_mu <- dm.jacobian(x, mu=as.numeric(mu1[n.step - 1]), 
                        alpha=alpha1[[n.step]], param = "mu")
    H_mu <- dm.hessian(x, mu=as.numeric(mu1[n.step - 1]), 
                       alpha=alpha1[[n.step]], param = "mu")
    
    ## Prevent NA or zero Hessian for mu
    if (is.na(H_mu) || H_mu == 0){
      mu1[n.step] <- mu1[n.step - 1]
    } else {
      W <- J_mu/H_mu
      step.size.check <- 10^(-(1:6))[
        sapply(1:6, \(x) as.numeric(mu1[n.step -1]) - (10^(-x))*W) >= 0
      ]
      step.size   <- if (length(step.size.check) == 0) 0 else step.size.check[1]
      mu1[n.step] <- mu1[n.step - 1] - step.size * W 
    }
    
    DM.lkhd[n.step] <- dm.lkhd(x, mu=as.numeric(mu1[n.step]), 
                               alpha=alpha1[[n.step]])
    eps <-  (DM.lkhd[n.step] - DM.lkhd[n.step - 1]) / sqrt(abs(DM.lkhd[n.step - 1]))
  }  
  
  # If last update is in the opposite direction, use penultimate update value             
  if(eps > 0 ){
    mu.est    <- mu1[n.step]
    alpha.est <- alpha1[n.step]
  } else {
    mu.est    <- mu1[n.step-1]
    alpha.est <- alpha1[n.step-1]
  }     
  
  return(list(mu=mu.est, alpha=alpha.est))
}
