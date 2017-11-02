
#' \code{S3} Class for functional datasets.
#' A class for univariate or multivariate functional dataset.
#' @param grid the grid (of length \code{T}) over which the functional dataset is defined.
#' @param data a vector, a matrix or a \code{list} of \code{d} vectors or matrices containing the functional data.
#'
#' @return The function returns a \code{S3} object of class \code{funData}, containing
#' the \code{grid} over which the functional dataset is defined and a matrix or a \code{list}
#' of \code{d} matrices containing the functional data
#'
#' @seealso \code{\link{simulate_KL}}
#'
#' @examples
#' # Define parameters
#' n <- 50
#' P <- 100
#' K <- 150
#'
#' # Grid of the functional dataset
#' t <- seq( 0, 1, length.out = P )
#'
#' # Define the means and the parameters to use in the simulation
#' # with the Karhunen - Loève expansion
#' m1 <- t^2 * ( 1 - t )
#' m2 <- t * ( 1 - t )^2
#'
#' lambda <- rep( 0, K )
#' theta <- matrix( 0, K, P )
#' for ( k in 1:K) {
#'   lambda[k] <- 1 / ( k + 1 )^2
#'   if ( k%%2 == 0 )
#'     theta[k, ] <- sqrt( 2 ) * sin( k * pi * t )
#'   else if ( k%%2 != 0 && k != 1 )
#'     theta[k, ] <- sqrt( 2 ) * cos( ( k - 1 ) * pi * t )
#'   else
#'     theta[k, ] <- rep( 1, P )
#' }
#' # Simulate the functional data
#' x1 <- simulate_KL( t, n, m1, rho = lambda, theta = theta )
#' x2 <- simulate_KL( t, n, m2, rho = lambda, theta = theta )
#'
#' FD <- funData( t, list( x1, x2 ) )

#' @export
#'
funData = function( grid, data )
{
  if (!is.list(data)) {
    return( structure( list( grid = grid, data = list(data) ),
                       class = c( 'funData' ) ) )
  }
  else {
    return( structure( list( grid = grid, data = data ),
                       class = c( 'funData' ) ) )
  }
}


#' A method to plot \code{funData} objects
#'
#' This function performs the plot of a functional univariate dataset stored in
#' an object of class \code{funData}. It is able to accept all the usual
#' customisable graphical parameters, otherwise it will use the default ones.
#'
#' @param x the univariate functional dataset in form of \code{funData} object.
#' @param ... additional graphical parameters to be used in plotting functions
#'
#' @seealso \code{\link{funData}}
#'
#' @examples
#'
#' N = 20
#'
#' @export
#'
#' @examples
#' # Define parameters
#' n <- 50
#' P <- 100
#' K <- 150
#'
#' # Grid of the functional dataset
#' t <- seq( 0, 1, length.out = P )
#'
#' # Define the means and the parameters to use in the simulation
#' # with the Karhunen - Loève expansion
#' m1 <- t^2 * ( 1 - t )
#' m2 <- t * ( 1 - t )^2
#'
#' lambda <- rep( 0, K )
#' theta <- matrix( 0, K, P )
#' for ( k in 1:K) {
#'   lambda[k] <- 1 / ( k + 1 )^2
#'   if ( k%%2 == 0 )
#'     theta[k, ] <- sqrt( 2 ) * sin( k * pi * t )
#'   else if ( k%%2 != 0 && k != 1 )
#'     theta[k, ] <- sqrt( 2 ) * cos( ( k - 1 ) * pi * t )
#'   else
#'     theta[k, ] <- rep( 1, P )
#' }
#' # Simulate the functional data
#' x1 <- simulate_KL( t, n, m1, rho = lambda, theta = theta )
#' x2 <- simulate_KL( t, n, m2, rho = lambda, theta = theta )
#'
#' FD <- funData( t, list( as.matrix( x1$data[[1]] ), as.matrix( x2$data[[1]] ) ) )
#'
#' plot(FD)
#'
plot.funData = function( x, ... ) {
  R <- length( x$data )
  par( mfrow = c( floor( sqrt( R ) ), ceiling( R / floor( sqrt( R ) ) ) ) )
  for ( r in 1:R ) {
    matplot( x$grid, t( x$data[[r]] ), type = "l", ... )
  }
}
