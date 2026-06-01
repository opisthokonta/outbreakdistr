


# Moment generating functions for infectious period distributions.
mgf <- function(t, params, model = "constant") {
  if (model == "constant") {
    # From Addy et al 1991 p 6.
    exp(-t * params)
  } else if (model == "exponential") {
    -params / (-params - t)
  } else if (model == "gamma") {
    (1 - (t / -params[2]))^(-params[1])
  }
}



validate_ip_model_params <- function(model, params) {
  if (!is.numeric(params)) {
    stop("Infectious period parameters must be numeric.")
  }

  if (!model %in% c("constant", "exponential", "gamma")) {
    stop('ip_model must be "constant", "exponential", or "gamma".')
  }

  if (model == "constant" | model == "exponential") {
    if (length(params) != 1) {
      stop(sprintf(
        "Number of parameters for infectious period must be 1 for %s", model
      ))
    }
    if (params <= 0) {
      stop("Infectious period parameter must be greater than 0.")
    }
  } else if (model == "gamma") {
    if (length(params) != 2) {
      stop(sprintf(
        "Number of parameters for infectious period must be 2 for %s", model
      ))
    }
  }
}


# Return a matrix where each row is one possible outcome (number infected per
# group). Number of columns equals number of groups.
make_multitype_state_table <- function(x) {
  stopifnot(all(x >= 0))
  m <- length(x)
  temp_list <- vector(mode = "list")
  for (ii in 1:m) {
    temp_list[[ii]] <- 0:x[ii]
  }
  statemat <- as.matrix(expand.grid(temp_list))
  colnames(statemat) <- NULL
  return(statemat)
}



# Return row indices of statemat where every element is <= the j-th row.
subset_statemat_idx <- function(statemat, j) {
  rows_to_keep <- c(rep(TRUE, j), rep(FALSE, nrow(statemat) - j))
  for (cc in 1:ncol(statemat)) {
    rows_to_keep <- rows_to_keep & (statemat[, cc] <= statemat[j, cc])
  }
  which(rows_to_keep)
}



# Compute the final size probability distribution (up to xmax infections).
# Based on Theorem 3.12 in Britton & Pardoux (2019).
fsdistr_internal <- function(s0, i0, beta, ip_model, ip_params, xmax = Inf) {
  stopifnot(
    s0 > 0,
    length(s0) == 1,
    i0 >= 1,
    length(i0) == 1,
    length(xmax) == 1,
    xmax >= 0,
    ip_model %in% c("constant", "exponential", "gamma")
  )

  xmax <- min(s0, xmax)
  N <- i0 + s0
  probs <- numeric(xmax + 1)

  for (kk in 0:xmax) {
    idx <- kk + 1

    if (kk == 0) {
      probs[idx] <- mgf(t = beta, params = ip_params, model = ip_model)^(i0)
      next
    }

    mgf_res <- mgf(
      t = ((N - i0 - kk) * (beta / s0)),
      params = ip_params,
      model = ip_model
    )

    term1 <- choose(N - i0, kk) * (mgf_res^(kk + i0))
    ii <- 0:(kk - 1)
    term2 <- sum(choose(N - i0 - ii, kk - ii) * (mgf_res^(kk - ii)) * probs[ii + 1])
    probs[idx] <- term1 - term2
  }

  if ((xmax == s0) & (abs(sum(probs) - 1) > 0.0001)) {
    warning("sum not 1")
  }
  if (any(probs < 0)) {
    warning("some probability negative")
  }

  return(probs)
}


#' Final size distribution for a single-type stochastic SIR epidemic
#'
#' Computes the exact probability distribution of the number of susceptibles
#' infected during an outbreak, for a closed homogeneously mixing population.
#' Based on Theorem 3.12 in Britton & Pardoux (2019).
#'
#' @param s0 Integer. Initial number of susceptibles (must be >= 1).
#' @param i0 Integer. Initial number of infectives (must be >= 1).
#' @param beta Numeric. Transmission rate parameter.
#' @param ip_model Character. Infectious period distribution: `"constant"`,
#'   `"exponential"`, or `"gamma"`.
#' @param ip_params Numeric vector. Parameters for the infectious period
#'   distribution. One value for `"constant"` or `"exponential"` (the rate or
#'   duration); two values `(shape, rate)` for `"gamma"`.
#'
#' @return A numeric vector of length `s0 + 1`. Element `k + 1` gives the
#'   probability that exactly `k` susceptibles are infected (k = 0, 1, ...,
#'   s0).
#'
#' @references
#' Britton, T., & Pardoux, E. (Eds.). (2019). *Stochastic Epidemic Models
#' with Inference*. Springer.
#'
#' @examples
#' # Probability distribution with exponential infectious period
#' fsdistr(s0 = 10, i0 = 1, beta = 1.5)
#'
#' # Constant infectious period of length 0.9
#' fsdistr(s0 = 10, i0 = 1, beta = 1.4, ip_model = "constant", ip_params = 1 / 0.9)
#'
#' @export
fsdistr <- function(s0, i0, beta, ip_model = "exponential", ip_params = 1) {
  validate_ip_model_params(model = ip_model, params = ip_params)
  stopifnot(length(ip_model) == 1)
  fsdistr_internal(
    s0 = s0, i0 = i0, beta = beta,
    ip_model = ip_model, ip_params = ip_params,
    xmax = Inf
  )
}


#' Final size distribution for a multi-type stochastic SIR epidemic
#'
#' Computes the exact joint probability distribution of the number infected
#' in each group during an outbreak, for a population with multiple
#' susceptible groups and heterogeneous mixing. Based on equation 3.4 in
#' Ball (1986).
#'
#' @param s0 Integer vector of length `m`. Initial susceptible counts per
#'   group.
#' @param i0 Integer vector of length `m`. Initial infective counts per
#'   group. At least one element must be >= 1.
#' @param beta Numeric matrix of dimension `m x m`. Entry `beta[i, j]` is
#'   the rate at which an infective in group `j` contacts susceptibles in
#'   group `i`.
#' @param ip_model Character vector of length `m` (or scalar, recycled).
#'   Infectious period model per group: `"constant"`, `"exponential"`, or
#'   `"gamma"`.
#' @param ip_params List of length `m` (or length 1, recycled). Each element
#'   is the parameter(s) for the infectious period distribution of the
#'   corresponding group.
#'
#' @return A named numeric vector of length `prod(s0 + 1)`. Each element
#'   gives the probability of the corresponding joint outcome, with names
#'   encoding the outcome as `"k1-k2-...-km"`.
#'
#' @references
#' Ball, F. (1986). A unified approach to the distribution of total size and
#' total area under the trajectory of infectives in epidemic models.
#' *Advances in Applied Probability*, 18(2), 289–310.
#'
#' @examples
#' # Two-group model
#' beta2 <- matrix(c(1.5, 0.5, 0.5, 1.5), nrow = 2)
#' fsdistr_mt(s0 = c(5, 5), i0 = c(1, 0), beta = beta2)
#'
#' # Single group (should match fsdistr())
#' fsdistr_mt(s0 = 6, i0 = 1, beta = matrix(1.4), ip_model = "constant",
#'            ip_params = list(1 / 0.9))
#'
#' @export
fsdistr_mt <- function(s0, i0, beta, ip_model = "exponential",
                       ip_params = list(1)) {
  m <- length(s0)

  if (length(ip_model) == 1) {
    ip_model <- rep(ip_model, m)
  }

  stopifnot(is.list(ip_params))
  if (length(ip_params) == 1) {
    ip_params <- rep(ip_params, m)
  }

  stopifnot(
    is.numeric(s0),
    all(s0 > 0),
    is.numeric(i0),
    all(i0 >= 0),
    any(i0 >= 1),
    length(i0) == m,
    is.matrix(beta),
    nrow(beta) == m,
    ncol(beta) == m,
    length(ip_model) == m,
    length(ip_params) == m
  )

  for (ii in 1:m) {
    validate_ip_model_params(model = ip_model[ii], params = ip_params[[ii]])
  }

  nstates <- prod(s0 + 1)
  statemat <- make_multitype_state_table(s0)
  statelabs <- apply(statemat, MARGIN = 1, FUN = \(x) paste(x, collapse = "-"))

  stopifnot(nrow(statemat) == nstates)

  probs <- numeric(nstates)

  for (jj in 1:nstates) {
    cur_state <- statemat[jj, ]
    yy <- prod(choose(s0, cur_state))

    statemat_idx <- subset_statemat_idx(statemat, j = jj)
    statematj <- statemat[statemat_idx, , drop = FALSE]

    xxj <- rep(0.0, jj)

    for (oo in 1:nrow(statematj)) {
      bcoef <- prod(choose(s0 - statematj[oo, ], cur_state - statematj[oo, ]))

      prod_tmp <- 1
      for (ii in 1:m) {
        mgf_eval_at <- sum((s0 - cur_state) * (beta[, ii] / s0[ii]))
        mgf_res_tmp <- mgf(
          t = mgf_eval_at,
          params = ip_params[[ii]],
          model = ip_model[ii]
        )^(statematj[oo, ii] + i0[ii])
        prod_tmp <- prod_tmp * mgf_res_tmp
      }

      xxj[statemat_idx[oo]] <- bcoef / prod_tmp
    }

    jidx <- 1:max(1, jj - 1)
    probs[jj] <- (yy - sum(probs[jidx] * xxj[jidx])) / xxj[jj]
  }

  names(probs) <- statelabs

  if (abs(sum(probs) - 1) > 0.0001) {
    warning("sum not 1")
  }
  if (any(probs < 0)) {
    warning("some probability negative")
  }

  return(probs)
}
