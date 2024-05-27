ortho.randproj <- function(nrow, ncol, method = "norm", seed = NULL){
    #' Orthogonal random matrix generator
    #'
    #' Generate random matrices of given dimension which have full row and are orthogonal
    #' @import mvtnorm stats
    NULL
    #' @export

    #' @param nrow Number of rows in the random matrix to be generated
    #' @param ncol Number of colums
    #' @param method The method to be used for generating elements in the matrix. If \code{method = "norm"} (Default), elements are generated from standard normal distribution. If \code{method = "achlioptas"}, the elements are selected from the set \eqn{\{-1, 0, 1\}} with probabilities \eqn{\{1/6, 2/3, 1/6\}} respectively.
    #' @param seed Set the seed to replicate the random matrix generated. Default is \code{NULL} and no seed is used.
    #' @return An orthogonal matrix of with \code{nrow} rows and \code{ncol} columns.
    
    #' @examples
    #' a = ortho.randproj(nrow = 5, ncol = 3)

    if (!is.null(seed)){
        set.seed(seed)
    }

    ## NORMAL random sample projection matrix
    if (method == "norm"){
        R <- matrix(rnorm(nrow*ncol), nrow = nrow, ncol = ncol)
    }

    ## Achlioptas projection matrix
    if (method == "achlioptas"){
        R <- matrix(sample( rep(c(-1,0,1)*sqrt(3), c(1,4,1)),  size = nrow*ncol, replace = TRUE), nrow = nrow, ncol = ncol)
    }

    ## ORTHOGONALIZE USING SVD
    W <- eigen(R%*%t(R))
    Rstar <- (W$vectors%*%( diag(1/sqrt(W$values))%*%t(W$vectors) ))%*%R

    return(Rstar)
}
