
# Enumerate all (S, I) states where S + I <= N.
make_markov_state_table <- function(N) {
  N_states <- sum(1:(N + 1))
  state_table <- matrix(nrow = N_states, ncol = 2)

  idx_start <- 1

  for (nn in N:0) {
    all_possible_s <- nn:0
    all_possible_i <- nn - all_possible_s

    idx_end <- idx_start + length(all_possible_s) - 1
    state_table[idx_start:idx_end, 1] <- all_possible_s
    state_table[idx_start:idx_end, 2] <- all_possible_i

    idx_start <- idx_end + 1
  }

  return(state_table)
}


# Look up the row index of a (S, I) state in the state table
get_state_table_idx <- function(sm, state) {
  stopifnot(ncol(sm) == 2, length(state) == 2)
  idx <- which(sm[, 1] == state[1] & sm[, 2] == state[2])
  stopifnot(length(idx) == 1)
  return(idx)
}

# Construct the transition matrix for the Epidemic Markov Model.
#
# statemat: matrix with one row per state, columns (S, I).
# s0:       initial number of susceptibles.
# beta:     infection rate.
# gamma:    recovery rate.
# beta0 :   seeding rate
# sparse:   if TRUE, return a Matrix::sparseMatrix; otherwise a dense matrix.
# elements: what kind of transition matrix to return, either "probabilities" or "rates".
#
# To obtain the canonical form, pass statemat with absorbing states first.
make_transition_matrix <- function(statemat, s0, beta, gamma, beta0 = 0, sparse = FALSE, elements = 'probabilities', model = 'sir') {

  if (model == 'sis' & elements == 'probabilities'){
    stop('Model SIS sohuld only be used with elements =  "rates".')
  }


  state_labels <- apply(statemat, MARGIN = 1, FUN = \(x) paste0(x, collapse = "-"))

  infection_rate <- (beta * statemat[, 1] * (statemat[, 2] / s0)) + (beta0 * statemat[, 1])
  recovery_rate <- gamma * statemat[, 2]

  if (elements == 'probabilities'){
    infection_probability <- infection_rate / (infection_rate + recovery_rate)
    recovery_probability <- recovery_rate / (infection_rate + recovery_rate)

    element_infection <- infection_probability
    element_recovery <- recovery_probability
    absorbing_state_element <- 1

  } else if (elements == 'rates'){
    element_infection <- infection_rate
    element_recovery <- recovery_rate
    absorbing_state_element <- 0
  }


  if (sparse) {

    if (model == 'sir'){
      next_state_inf_labels <- paste(statemat[, 1] - 1, statemat[, 2] + 1, sep = "-")
      next_state_rec_labels <- paste(statemat[, 1], statemat[, 2] - 1, sep = "-")
    } else if (model == 'sis'){
      next_state_inf_labels <- paste(statemat[, 1] - 1, statemat[, 2] + 1, sep = "-")
      next_state_rec_labels <- paste(statemat[, 1] + 1, statemat[, 2] - 1, sep = "-")
    }


    next_state_inf_idx <- match(next_state_inf_labels, state_labels)
    next_state_rec_idx <- match(next_state_rec_labels, state_labels)

    absorbing_state <- statemat[, 2] == 0
    absorbing_state_idx <- which(absorbing_state)

    # Transient states where both infection and recovery are possible.
    transient_with_s <- (!absorbing_state) & (statemat[, 1] >= 1)
    # Transient states where only recovery is possible.
    transient_without_s <- (!absorbing_state) & (statemat[, 1] == 0)

    transient_with_s_idx <- which(transient_with_s)
    transient_without_s_idx <- which(transient_without_s)

    rows <- c(
      absorbing_state_idx,
      transient_without_s_idx,
      transient_with_s_idx,
      transient_with_s_idx
    )

    cols <- c(
      absorbing_state_idx,
      next_state_rec_idx[transient_without_s],
      next_state_rec_idx[transient_with_s],
      next_state_inf_idx[transient_with_s]
    )

    vals <- c(
      rep(absorbing_state_element, length(absorbing_state_idx)),
      element_recovery[transient_without_s],
      element_recovery[transient_with_s],
      element_infection[transient_with_s]
    )

    transition_matrix <- Matrix::sparseMatrix(
      i = rows, j = cols, x = vals,
      dimnames = list(state_labels, state_labels)
    )
  } else { # Dense.

    transition_matrix <- matrix(
      0,
      ncol = nrow(statemat),
      nrow = nrow(statemat),
      dimnames = list(state_labels, state_labels)
    )

    if (model == 'sir'){

      for (ii in 1:nrow(statemat)) {
        cur_state <- statemat[ii, ]
        cur_state_idx <- ii

        if (cur_state[2] == 0) {
          # Absorbing state: stays in place.

          next_state_idx <- get_state_table_idx(sm = statemat, state = cur_state)
          transition_matrix[cur_state_idx, next_state_idx] <- absorbing_state_element # TODO: what if beta0 != 0??

        } else if (cur_state[2] > 0 & cur_state[1] > 0) {
          # Infection or recovery both possible.
          next_state_inf <- c(cur_state[1] - 1, cur_state[2] + 1)
          next_state_idx <- get_state_table_idx(sm = statemat, state = next_state_inf)
          transition_matrix[cur_state_idx, next_state_idx] <- element_infection[ii]

          next_state_rec <- c(cur_state[1], cur_state[2] - 1)
          next_state_idx <- get_state_table_idx(sm = statemat, state = next_state_rec)

          transition_matrix[cur_state_idx, next_state_idx] <- element_recovery[ii]

        } else if (cur_state[2] > 0 & cur_state[1] == 0) {
          # Only recovery possible (no susceptibles left).

          next_state_rec <- c(cur_state[1], cur_state[2] - 1)
          next_state_idx <- get_state_table_idx(sm = statemat, state = next_state_rec)
          transition_matrix[cur_state_idx, next_state_idx] <- element_recovery[ii]
        }
      }
    } else if (model == 'sis'){

      for (ii in 1:nrow(statemat)) {
        cur_state <- statemat[ii, ]
        cur_state_idx <- ii

        if (cur_state[2] > 0){
          # If there are any infectious
          # they can recover back to be susceptible

          next_state <- c(cur_state[1] + 1,  cur_state[2] - 1)
          next_state_idx <- get_state_table_idx(sm = statemat, state = next_state)

          transition_matrix[cur_state_idx, next_state_idx] <- recovery_rate[ii]

        }

        if (cur_state[2] > 0 & cur_state[1] > 0){
          # When there any infectious AND any susceptible, then a susceptible can be infected.

          next_state <- c(cur_state[1] - 1,  cur_state[2] + 1)
          next_state_idx <- get_state_table_idx(sm = statemat, state = next_state)

          transition_matrix[cur_state_idx, next_state_idx] <- infection_rate[ii]
        }
      }
    }
  }

  if (elements == 'probabilities'){
    row_sums <- if (sparse) Matrix::rowSums(transition_matrix) else rowSums(transition_matrix)
    stopifnot(all(abs(row_sums - 1) <= 0.00001))
  }

  return(transition_matrix)
}





#' Epidemic Markov Models
#'
#' Exact computation for Epidemic Markov Models. Use the `emmct` for continuous time
#' modeling of a SIR epidemic, use `emmct_sis` for epidemics with a SIS dynamic.
#' `emmdt` is the discrete time variant of the SIR epidemic, which is mostly useful
#' for computing the final size distribution. The models are described in detail in
#' the Epidemic Markov Models vignette.
#'
#'
#' @param s0 Positive integer. Initial number of susceptibles.
#' @param i0 Positive integer. Initial number of infectives (must be >= 1).
#' @param beta Non-negative numeric. Infection rate parameter. The per-step
#'   infection probability is proportional to `beta * S * I / s0`.
#' @param gamma Non-negative numeric. Recovery rate parameter. The per-step
#'   recovery probability is proportional to `gamma * I`.
#' @param time Positive numeric. The time to evaluate the the model, assuming an initial state at time 0 defined by `s0` and `i0`.
#' @param beta0 Non-negative numeric. Rate of infection from outside the population.
#' @param sparse Logical. If `TRUE` (default), use a sparse matrix
#'   representation via the Matrix package; otherwise use a dense matrix.
#'   Sparse matrices are faster for large populations.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{`transition_matrix`}{The full one-step transition matrix.}
#'   \item{`qmat`}{Sub-matrix of transitions among transient states
#'     (canonical form).}
#'   \item{`rmat`}{Sub-matrix of transitions from transient to absorbing
#'     states (canonical form).}
#'   \item{`fmat`}{Fundamental matrix `(I - Q)^{-1}`.}
#'   \item{`solution_mat`}{Absorption probability matrix: entry `[i, j]` is
#'     the probability of being absorbed in state `j` when starting in
#'     transient state `i`.}
#'   \item{`fs_distr`}{Numeric vector of length `s0 + 1`. Element `k + 1`
#'     gives the probability that exactly `k` susceptibles are infected
#'     (k = 0, 1, ..., s0).}
#' }
#'
#' @examples
#' res <- emmdt(s0 = 10, i0 = 1, beta = 1.5, gamma = 1.0)
#' res$fs_distr
#'
#' # Dense matrix representation
#' res_dense <- emmdt(s0 = 10, i0 = 1, beta = 1.5, gamma = 1.0, sparse = FALSE)
#'
#' @export
emmdt <- function(s0, i0, beta, gamma, sparse = TRUE) {

  stopifnot(
    length(s0) == 1,
    length(i0) == 1,
    s0 >= 1,
    i0 >= 1,
    beta >= 0,
    gamma >= 0,
    is.logical(sparse)
  )

  N <- s0 + i0
  state_table <- make_markov_state_table(N)

  transition_matrix <- make_transition_matrix(
    state_table,
    s0 = s0,
    beta = beta,
    gamma = gamma,
    sparse = sparse,
    elements = 'probabilities'
  )

  absorbing_states <- state_table[, 2] == 0
  absorbing_states_idx <- which(absorbing_states)
  transient_states_idx <- which(!absorbing_states)

  qmat <- transition_matrix[transient_states_idx, transient_states_idx]
  rmat <- transition_matrix[transient_states_idx, absorbing_states_idx]

  # Compute the fundamental matrix
  if (sparse){
    fmat <- Matrix::solve(Matrix::Diagonal(ncol(qmat)) - qmat)
  } else {
    fmat <- solve(diag(ncol(qmat)) - qmat)
  }

  # Compute the solution matrix, i.e. probability of of absorption in state j while
  # starting in state i.
  solution_mat <- fmat %*% rmat

  if (sparse){
    sol_row_sums <- Matrix::rowSums(solution_mat)
  } else {
    sol_row_sums <- rowSums(solution_mat)
  }


  stopifnot(
    all(abs(sol_row_sums - 1) <= 0.0001),
    all(solution_mat >= 0)
  )

  # Final size distribution
  i_final <- s0 - state_table[absorbing_states_idx, 1]

  # Verify absorbing states are ordered as expected.
  stopifnot(all(i_final == (-i0:s0)))

  fs_distr <- solution_mat[paste0(c(s0, i0), collapse = "-"), i_final >= 0]
  names(fs_distr) <- NULL

  stopifnot(all(abs(sum(fs_distr) - 1) <= 0.0001),
            all(fs_distr >= 0),
            length(fs_distr) == s0+1)


  res <- list(transition_matrix = transition_matrix,
              qmat = qmat,
              rmat = rmat,
              fmat = fmat,
              solution_mat = solution_mat,
              fs_distr = fs_distr
              )

  return(res)

}




#' @rdname emmdt
#' @export
emmct <- function(s0, i0, beta, gamma, time, beta0 = 0, sparse = TRUE, expm_method = NULL){

  stopifnot(
    length(s0) == 1,
    length(i0) == 1,
    s0 >= 1,
    i0 >= 1,
    beta >= 0,
    gamma >= 0,
    is.logical(sparse)
  )

  if (is.null(expm_method)){
    if (sparse){
      expm_method <- "Higham08"
    } else {
      expm_method <- "Higham08.b"
    }
  } else {
    stopifnot(length(expm_method) == 1,
              is.character(expm_method))
  }

  N <- s0 + i0
  state_table <- make_markov_state_table(N)

  qmatrix <- make_transition_matrix(statemat = state_table,
                                    s0 = s0,
                                    beta = beta,
                                    gamma = gamma,
                                    beta0 = beta0,
                                    sparse = sparse,
                                    elements = 'rates')


  if (sparse){
    diag(qmatrix) <- -Matrix::rowSums(qmatrix)
  } else {
    diag(qmatrix) <- -rowSums(qmatrix)
  }


  # Compute the transition matrix.
  transition_matrix <- expm::expm(qmatrix * time, method = expm_method)

  # Small check.
  if (sparse){
    stopifnot(all(abs(Matrix::rowSums(transition_matrix) - 1) <= 0.0001))
  } else {
    stopifnot(all(abs(rowSums(transition_matrix) - 1) <= 0.0001))
  }


  # Identify the row corresponding to the initial state.
  initial_state_idx <- which(state_table[,1] == s0 & state_table[,2] == i0)

  state_prob <- transition_matrix[initial_state_idx,]

  # The cumulative number of infected, for each state.
  i_t <- s0 - state_table[,1]

  # Compute cumulative incidence probability distribution
  # (ie  number of susceptibles who has become infected).
  fs_distr <- numeric(s0+1)
  for (ii in 0:(length(fs_distr)-1)){
    fs_distr[ii+1] <- sum(state_prob[i_t == ii])
  }

  stopifnot(all(abs(sum(fs_distr) - 1) <= 0.0001),
            all(fs_distr >= 0),
            length(fs_distr) == s0+1)


  # Prepare output
  res <- list(qmatrix = qmatrix,
              transition_matrix = transition_matrix,
              time = time,
              state_prob = state_prob,
              fs_distr = fs_distr)

  return(res)

}




#' @rdname emmdt
#' @export
emmct_sis <- function(s0, i0, gamma, beta, time, beta0 = 0, sparse = TRUE, expm_method = NULL){

  if (is.null(expm_method)){
    if (sparse){
      expm_method <- "Higham08"
    } else {
      expm_method <- "Higham08.b"
    }
  } else {
    stopifnot(length(expm_method) == 1,
              is.character(expm_method))
  }

  N <- s0 + i0

  #  Make the SIS state table.
  state_table <- matrix(nrow = N+1, ncol = 2)

  state_table[,1] <- N:0
  state_table[,2] <- 0:N

  qmatrix <- make_transition_matrix(statemat = state_table,
                                    s0 = s0,
                                    beta = beta,
                                    gamma = gamma,
                                    beta0 = beta0,
                                    sparse = sparse,
                                    elements = 'rates',
                                    model = 'sis')


  if (sparse){
    diag(qmatrix) <- -Matrix::rowSums(qmatrix)
  } else {
    diag(qmatrix) <- -rowSums(qmatrix)
  }

  # Compute the transition matrix.
  transition_matrix <- expm::expm(qmatrix * time, method = expm_method)

  # Small check.
  if (sparse){
    stopifnot(all(abs(Matrix::rowSums(transition_matrix) - 1) <= 0.0001))
  } else {
    stopifnot(all(abs(rowSums(transition_matrix) - 1) <= 0.0001))
  }

  # Identify the row corresponding to the initial state.
  initial_state_idx <- which(state_table[,1] == s0 & state_table[,2] == i0)

  state_prob <- transition_matrix[initial_state_idx,]

  # Prepare output
  res <- list(qmatrix = qmatrix,
              transition_matrix = transition_matrix,
              time = time,
              state_prob = state_prob)

  return(res)

}





