random.proj <- function(nrow, ncol, seed=NULL){ 
  #' @title Random matrix generator
  #' 
  #' @description This function generates a random matrix whose elements are zeros and ones. The matrix is of full row rank.
  #' @param nrow Number of rows in the random matrix to be generated
  #' @param ncol Number of columns
  #' @param seed Seed for random number generation. Default is NULL, where a random seed is used. However, the seed can be specified for reproducibility.
  #' @import stats
  #' @export
  #' @return A matrix of specified dimension with zero-one switch based random projections. The matrix is verified to have full row rank.
  
  R <- matrix(0, nrow = nrow, ncol = ncol)
  
  repeat{
    for (j in 1:ncol){
      R[,j]<- rmultinom(1, 1, prob = rep(1/nrow, nrow))
    }
    if(qr(R)$rank == nrow){   
      break
    }} 
  return(R)
}

