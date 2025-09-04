#' Function to calculate sample size corresponding to desired power
#'
#' @param power Numeric; desired level of power
#' @param ... Additional arguments to pass to swCRTdesign::swPwr
#' @return A list containing keys 'n' and 'power' corresponding to the sample
#'     size (per cluster-period) and exact power
#'     Author: Avi Kenny    7/5/24
swSiz <- function(power, alpha=0.05, ...) {
  
  if (power <= alpha | power >= 1) stop("power must be greater than alpha (default = 0.05) and less than 1")

    power_n <- function(n) {
    o_args <- list(...)
    o_args$alpha <- alpha
    o_args$n <- n
    o_args$retDATA <- FALSE
    o_args$silent <- TRUE
    suppressMessages({ do.call(swCRTdesign::swPwr, o_args) })
  }  
  
  # Set initial sample size range and corresponding power
  n_lo <- 1
  n_up <- 2
  power_lo <- power_n(n_lo)
  power_up <- power_n(n_up)
  if (power_lo>power) { stop("Power achieved with n=1.") }
  
  # Get initial sample size range that contains 90% power
  while(power_up<power) {
    n_lo <- n_up
    power_lo <- power_up
    n_up <- 2*n_up
    power_up <- power_n(n_up)
  }
  
  # Bisect interval until 90% range is found
  while(n_up-n_lo>1) {
    n_mid <- round(mean(c(n_lo,n_up)))
    power_mid <- power_n(n_mid)
    if (power_mid<power) {
      n_lo <- n_mid
      power_lo <- power_mid
    } else {
      n_up <- n_mid
      power_up <- power_mid
    }
  }
  
  if (abs(power-power_lo)<abs(power-power_up)) {
    return(list(n=n_lo, power=power_lo))
  } else {
    return(list(n=n_up, power=power_up))
  }
  
}


