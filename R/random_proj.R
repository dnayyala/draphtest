
random.proj <- function(nrow, ncol, seed=NULL){ 
#' Random matrix generator for Dirichlet distribution
#' 
#' @param nrow Number of rows in the random matrix to be generated
#' @param ncol Number of columns
#' @param method The method to be used for generating elements in the matrix. If   


  R <- matrix(NA, k, p)
  
  repeat{
    for (j in 1:p){
      R[,j]<- rmultinom(1, 1, prob = rep(1/k, k))
    }
    if(qr(R)$rank == k){   
      break
    }} 
  return(R)
}

