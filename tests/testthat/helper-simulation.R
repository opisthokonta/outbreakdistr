# Simulation helpers for testing fsdistr() and fsdistr_mt() against Monte Carlo results.
# This file is automatically sourced by testthat before running tests.


# Removes Inf and -Inf elements in a vector.
remove_inf <- function(x) {
  x[!x %in% c(Inf, -Inf)]
}


# Counts occurrences of each group index 1:m in vector x.
count_i <- function(x, m) {
  res <- integer(m)
  for (ii in seq_len(m)) {
    res[ii] <- sum(x == ii)
  }
  res
}


# rexp wrapper that returns NA (not NaN, no warning) when rate = 0.
rexp0 <- function(n, rate) {
  if (any(rate == 0)) {
    is0 <- rate == 0
    res <- numeric(n)
    res[is0]  <- NA_real_
    res[!is0] <- rexp(n = sum(!is0), rate = rate[!is0])
  } else {
    res <- rexp(n, rate)
  }
  res
}


# Build the expansion-grid state matrix (matches make_statemat / make_multitype_state_table).
make_statemat_sim <- function(x) {
  stopifnot(all(x >= 0))
  temp_list <- lapply(x, function(xi) 0:xi)
  statemat  <- as.matrix(expand.grid(temp_list))
  colnames(statemat) <- NULL
  statemat
}


# Draw one infectious-period realisation.
simulate_infectious_period <- function(model, params) {
  if (model == "constant") {
    params
  } else if (model == "exponential") {
    rexp(1, rate = params)
  } else if (model == "gamma") {
    rgamma(1, shape = params[1], rate = params[2])
  }
}


# Simulate one SIR epidemic (possibly multi-type).
# Returns the number of new infections per group (integer vector of length m).
run_multitype_simulation <- function(i0, s0, beta,
                                     ip_model  = "exponential",
                                     ip_params = list(1)) {
  m <- length(s0)

  if (length(ip_model) == 1) ip_model <- rep(ip_model, m)

  stopifnot(is.list(ip_params))
  if (length(ip_params) == 1) ip_params <- rep(ip_params, m)

  stopifnot(
    is.numeric(s0), all(s0 > 0),
    is.numeric(i0), all(i0 >= 0), any(i0 >= 1),
    length(i0) == m,
    is.matrix(beta), nrow(beta) == m, ncol(beta) == m,
    length(ip_model) == m,
    length(ip_params) == m
  )

  N <- sum(s0 + i0)
  individual_group_idx <- rep(seq_len(m), times = s0 + i0)
  infection_times  <- rep(Inf, N)
  recovery_times   <- rep(Inf, N)

  # Seed initial infectives
  for (mm in seq_len(m)) {
    if (i0[mm] > 0) {
      for (ii in seq_len(i0[mm])) {
        idx <- min(which(individual_group_idx == mm)) + ii - 1
        infection_times[idx] <- 0
        recovery_times[idx]  <- simulate_infectious_period(ip_model[mm], ip_params[[mm]])
      }
    }
  }

  cumulative_i <- i0
  cur_s    <- s0
  cur_time <- 0

  while (any(cur_s > 0)) {

    cur_i_idx <- which(infection_times <= cur_time & recovery_times > cur_time)
    if (length(cur_i_idx) == 0) break

    cur_i_group <- count_i(individual_group_idx[cur_i_idx], m = m)

    foi              <- as.numeric(cur_i_group %*% (beta / s0))
    infection_rate   <- foi * cur_s
    next_inf_times   <- cur_time + rexp0(n = m, rate = infection_rate)
    next_inf_time    <- min(next_inf_times, na.rm = TRUE)
    next_inf_group   <- which.min(next_inf_times)

    rt1          <- remove_inf(recovery_times)
    rt1          <- rt1[rt1 > cur_time]
    it1          <- remove_inf(infection_times)
    it1          <- it1[it1 > cur_time]
    next_event   <- min(rt1, it1)

    if (next_inf_time < next_event) {
      idx_to_update <- which(
        individual_group_idx == next_inf_group
      )[cumulative_i[next_inf_group] + 1]

      stopifnot(length(idx_to_update) == 1, idx_to_update <= N)

      infection_times[idx_to_update] <- next_inf_time
      recovery_times[idx_to_update]  <- next_inf_time +
        simulate_infectious_period(ip_model[next_inf_group], ip_params[[next_inf_group]])

      cur_s[next_inf_group]        <- cur_s[next_inf_group] - 1
      cumulative_i[next_inf_group] <- cumulative_i[next_inf_group] + 1
    }

    cur_time <- min(next_event, next_inf_time)
  }

  stopifnot(
    sum(cumulative_i) <= N,
    all(cur_s >= 0),
    sum(cumulative_i) == length(infection_times[infection_times < Inf])
  )

  cumulative_i - i0
}


# Simulate the single-type final-size distribution (returns probability vector
# of length s0 + 1, indexed 0:s0).
simulation_distr <- function(nsim, i0, s0, beta, ip_model, ip_params) {
  stopifnot(length(s0) == 1)

  counts <- integer(s0 + 1)
  for (ii in seq_len(nsim)) {
    res <- run_multitype_simulation(
      i0 = i0, s0 = s0, beta = beta,
      ip_model = ip_model, ip_params = ip_params
    )
    counts[res + 1] <- counts[res + 1] + 1L
  }
  counts / nsim
}



# Simulate the multi-type final-size distribution.
# Returns the statemat with an extra column of counts.
simulation_distr_mt <- function(nsim, i0, s0, beta, ip_model, ip_params) {
  m        <- length(s0)
  statemat <- cbind(make_statemat_sim(s0), 0L)
  statelabs <- apply(
    statemat[, seq_len(m), drop = FALSE], 1,
    function(x) paste(x, collapse = "-")
  )

  for (ii in seq_len(nsim)) {
    sim_res   <- run_multitype_simulation(
      i0 = i0, s0 = s0, beta = beta,
      ip_model = ip_model, ip_params = ip_params
    )
    state_idx <- match(paste(sim_res, collapse = "-"), statelabs)
    statemat[state_idx, m + 1] <- statemat[state_idx, m + 1] + 1L
  }

  stopifnot(sum(statemat[, m + 1]) == nsim)
  statemat
}
