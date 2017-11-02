#' Two-sample hypotesis tests
#'
#' Performs a  two sample hypotesis tests on two samples of functional data.
#' @param grid the grid (of length \code{T}) over which the functional dataset is defined.
#' @param x a matrix or a \code{list} of \code{d} matrices containing the first sample
#' @param y a matrix or a \code{list} of \code{d} matrices containing the second sample
#' @param conf.level confidence level of the test
#' @param stat_test the chosen test statistic to be used: \code{"L2"} for the classical L2-distance, \code{"L2_trunc"} for the truncated L2-distance, \code{"trunc"} for the truncated Mahalanobis semi-distance, \code{"mahalanobis"} for the generalized Mahalanobis distance
#' @param p a vector of positive numeric value containing the parameters of the regularizing function for the generalized Mahalanobis distance.
#' @keywords Inference
#' @return The function returns a list with the following components:
#'
#' \code{statistic} the value of the test statistic.
#'
#' \code{p.value} the p-value for the test.
#' @references
#'
#' Ghiglietti A., Ieva F., Paganoni A. M. (2017). Statistical inference for stochastic processes:
#' Two-sample hypothesis tests, \emph{Journal of Statistical Planning and Inference}, 180:49-68.
#'
#' Ghiglietti A., Paganoni A. M. (2017). Exact tests for the means of gaussian stochastic processes.
#' \emph{Statics & Probability Letters}, 131:102--107.
#'
#' @seealso \code{\link{funDist}}
#' @import
#' graphics
#' stats
#' @export
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
#'
#' s <- 0
#' for (k in 4:K) {
#'  s <- s + sqrt(lambda[k])*theta[k,]
#' }
#'
#' m2 <- m1 + s
#'
#' # Simulate the functional data
#' x1 <- simulate_KL( t, n, m1, rho = lambda, theta = theta )
#' x2 <- simulate_KL( t, n, m2, rho = lambda, theta = theta )
#'
#' output <- dp.test(t, x1$data, x2$data, 0.95, "L2", 1)

dp.test <- function( grid, x, y, conf.level = 0.95, stat_test, p=NULL ) {
  if ( !is.list( x ) ) {
    x <- list( x )
  }
  if ( !is.list( y ) ) {
    y <- list( y )
  }
  n1tot <- dim( x[[1]] )[1]
  n2tot <- dim( y[[1]] )[1]
  Ntot <- n1tot + n2tot
  n1 <- round( 0.75 * n1tot )
  n2 <- n1
  N <- n1 + n2
  P <- length( grid )
  R <- length( x )
  I <- grid[P]
  x_red <- list( )
  y_red <- list( )
  ind1 <- sample( 1:n1tot, n1, replace=F )
  ind2 <- sample( 1:n2tot, n2, replace=F )
  for ( r in 1:R ) {
    x_red[[r]] <- x[[r]][ind1, ]
    y_red[[r]] <- y[[r]][ind2, ]
  }

  xm <- list( )
  ym <- list( )
  for ( r in 1:R ) {
    xm[[r]] <- colMeans( x_red[[r]] )
    ym[[r]] <- colMeans( y_red[[r]] )
  }

  M <- min( N - 1, P - 1 )				# numero di autovalori stimati non nulli
  k_trunc <- 2
  xR <- NULL
  yR <- NULL
  for ( r in 1:R ) {
    xR <- cbind( xR, x_red[[r]] )
    yR <- cbind( yR, y_red[[r]] )
  }

  v1hat <- cov( xR )
  v2hat <- cov( yR )
  cn <- n1 / N
  vhat <- ( 1 - cn ) * v1hat + cn * v2hat
  Vhat <- eigen( vhat )$vectors * sqrt( P / I )	# autofunzioni di vhat
  dhat <- eigen( vhat )$values * I / P		# autovalori di vhat

  iter <- 100
  T0_1pop_vhat <- rep( 0, iter )	# le statistiche test con i diversi p
  qv <- rep( 0, length( p ) )			# i quantili con i diversi p

  ### per ogni p, trovo i quantili simulando la distribuzione della statistica test sotto H0, con gli autovalori stimati

  dist_H0_vhat <- rep( 0, 10^3 )		# Bootstrap sulle distanze sotto H0 per calcolare il quantile

  if ( stat_test == "mahalanobis" ) {
    for ( w in 1:10^3 ) {
      dist_H0_vhat[w] <- sum( rchisq( M, 1 ) * hhat( p, dhat[1:M] ) )
    }
  }
  #else if (stat_test == "L2_trunc") {
  #  for (w in 1:10^3) {
  #    dist_H0_vhat[w] <- sum(rchisq(k_trunc,1))
  #  }
  #}
  else if ( stat_test == "L2" ) {
    for ( w in 1:10^3 ) {
      dist_H0_vhat[w] <- sum( rchisq( M, 1 ) * dhat[1:M] )
    }
  }
  else if ( stat_test == "trunc" ) {
    for ( w in 1:10^3 ) {
      dist_H0_vhat[w] <- sum( rchisq( k_trunc, 1 ) * dhat[1:k_trunc] )
    }
  }

  qv <- as.numeric( quantile( dist_H0_vhat, conf.level ) )

  for ( j in 1:iter ) {
    x_red <- list( )
    y_red <- list( )
    ind1 <- sample( 1:n1tot, n1, replace = F )
    ind2 <- sample( 1:n2tot, n2, replace = F )
    for ( r in 1:R ) {
      x_red[[r]] <- x[[r]][ind1, ]
      y_red[[r]] <- y[[r]][ind2, ]
    }
    xm <- list( )
    ym <- list( )
    for ( r in 1:R ) {
      xm[[r]] <- colMeans( x_red[[r]] )
      ym[[r]] <- colMeans( y_red[[r]] )
    }

    xR <- NULL
    yR <- NULL
    for ( r in 1:R ) {
      xR <- cbind( xR, x_red[[r]] )
      yR <- cbind( yR, y_red[[r]] )
    }

    v1hat <- cov( xR )
    v2hat <- cov( yR )
    cn <- n1 / N
    vhat <- ( 1 - cn ) * v1hat + cn * v2hat
    Vhat <- eigen( vhat )$vectors# * sqrt( P / I )	# autofunzioni di vhat
    dhat <- eigen( vhat )$values #* I / P		# autovalori di vhat

    if ( stat_test == "mahalanobis" ) {
      T0_1pop_vhat[j] <- ( 1 / n1 + 1 / n2 )^( - 1 / 2 ) * funDist( grid, xm, ym, metric = "mahalanobis", lambda = dhat, phi = Vhat, p = p )
    }
    else if ( stat_test == "L2" ) {
      T0_1pop_vhat[j] <- ( 1 / n1 + 1 / n2 )^( - 1 / 2 ) * funDist( grid, xm, ym, metric = "L2" )
    }
    else if ( stat_test == "trunc" ) {
      T0_1pop_vhat[j] <- ( 1 / n1 + 1 / n2 )^( - 1 / 2 ) * funDist( grid, xm, ym, metric = "trunc", lambda = dhat, phi = Vhat )
    }
  }

  colMeans( T0_1pop_vhat > cbind( rep( 1, iter ) ) %*% rbind( qv ) )

  list ( T0 = T0_1pop_vhat,
         p = colMeans( T0_1pop_vhat > cbind( rep( 1, iter ) ) %*% rbind( qv ) )
  )
}
