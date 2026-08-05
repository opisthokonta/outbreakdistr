

# Internal functions to be used by get_transmission_probability.
# Infectious period functions to integrate.
# ie the expected value of exp(-(beta/s0)*x)
intf_exponential <- function(x, beta, s0, ip_params){
  exp(-(beta/s0) * x) * dexp(x = x, rate = ip_params)
}

intf_gamma <- function(x, beta, s0, ip_params){
  exp(-(beta/s0) * x) * dgamma(x = x, shape = ip_params[1], rate = ip_params[2])
}


#' Transmission probability for a stochastic SIR epidemic
#'
#' Computes the probability that a single infected individual infects a given
#' susceptible during their infectious period, in a stochastic SIR model.
#' Supports constant, exponential, and gamma infectious period distributions.
#'
#' @param beta Numeric vector. Transmission rate parameter(s).
#' @param s0 Numeric vector. Initial number of susceptibles. Must be the same
#'   length as `beta`.
#' @param ip_model Character. Infectious period distribution: `"constant"`,
#'   `"exponential"`, or `"gamma"`.
#' @param ip_params Numeric vector. Parameters for the infectious period
#'   distribution. One value for `"constant"` or `"exponential"` (the rate or
#'   duration); two values `(shape, rate)` for `"gamma"`.
#'
#' @return A numeric vector of the same length as `beta`, with each element
#'   giving the transmission probability for the corresponding `beta` and `s0`
#'   values. Values are in the interval \eqn{[0, 1]}.
#'
#' @examples
#' transmission_probability(beta = 1.1, s0 = 1, ip_model = "constant", ip_params = 1)
#' transmission_probability(beta = 1.1, s0 = 1, ip_model = "exponential", ip_params = 1)
#' transmission_probability(beta = 1.1, s0 = 1, ip_model = "gamma", ip_params = c(2, 2))
#'
#' @export
transmission_probability <- function(beta, s0, ip_model = 'constant', ip_params = 1){

  stopifnot(length(beta) == length(s0))

  res <- numeric(length(beta))

  for (ii in 1:length(beta)){
    if (ip_model == 'constant'){
      stopifnot(length(ip_params) == 1)
      res[ii] <- 1 - exp(-(beta[ii]/s0[ii])*ip_params)
    } else if (ip_model == 'exponential'){
      stopifnot(length(ip_params) == 1)
      upper_lim <- qexp(0.999, rate = ip_params)
      integrate_res <- integrate(f = intf_exponential,
                                 lower = 0, upper = upper_lim,
                                 beta = beta[ii], s0 = s0[ii],
                                 ip_params = ip_params)
      res[ii] <- 1 - integrate_res$value
    } else if (ip_model == 'gamma'){
      stopifnot(length(ip_params) == 2)
      upper_lim <- qgamma(0.999, shape = ip_params[1],  rate = ip_params[2])
      integrate_res <- integrate(f = intf_gamma,
                                 lower = 0, upper = upper_lim,
                                 beta = beta[ii], s0 = s0[ii],
                                 ip_params = ip_params)
      res[ii] <- 1 - integrate_res$value
    }
  }

  return(res)

}





