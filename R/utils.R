integral <- function ( grid, values ) {
  return( sum( ( values[ -1 ] + values[ -length( values ) ] ) * diff( grid ) / 2 ) )
}

#' Regularizing h function
#'
#' This function computes the values of the regularizing function for the generalized Mahalanobis distance.
#' @param p a positive numeric value containing the parameter of the regularizing function for the generalized Mahalanobis distance.
#' @param lambda a vector containing the eigenvalues of the covariance matrix of the functional data from which \code{"x"} and \code{"y"} are extracted.
#' @return The function returns a vector of the values of the regularizing function for the generalized Mahalanobis distance.
#' @export
#' @examples
#' # Define the parameters of the function
#' p <- 10^3
#' K <- 150
#' lambda <- rep( 0, K )
#' for ( k in 1:K ) {
#'   lambda[k] <- 1 / ( k + 1 )^2
#' }
#'
#' h <- hhat( p, lambda )

hhat <- function( p , lambda ) {  	# definisco la funzione h
  return( lambda / ( lambda + 1 / p ) )
}
