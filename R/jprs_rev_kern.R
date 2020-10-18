#' Single tour of the Jump Process Restore Sampler when the jump chain of
#' the underlying jump process has a reversible transition kernel.
#'
#' Assumes the holding rate is 1.
#' Assumes C * (normalising constant of target) = 1.
#'
#' @param dtarg_log Function to evaluate the log target density
#' @param dmu_log Function to evaluate the log regeneration distribution
#' @param rmu Function to simulate n times from the regeneration distribution
#' @param mh_kern Metropolis-Hastings kernel.
#' Input: current state. Output next state.
#' @return A matrix with first d rows (d the dimension of the target density)
#' the jump chain and last row the holding time for the jump process.
jprs_tour <- function(dtarg_log, dmu_log, rmu, mh_kern){

  # Regenerate
  x0 <- rmu(1)
  d <- length(x0)

  # Matrix to store samples
  jc_max_len <- 100
  x <- matrix(0, nrow=jc_max_len, ncol=(d+1))

  i <- 1
  x[1,1:d] <- x0
  regen_yet_to_occur <- 1

  while (regen_yet_to_occur){

    # Check x has enough rows:
    if ((i+1) > jc_max_len){
      # Double storage capacity if needed
      rbind(x, matrix(0, nrow=jc_max_len, ncol=(d+1)))
      jc_max_len <- jc_max_len + jc_max_len
    }

    # Time of local move
    t_local <- rexp(1, rate=1)

    # Regeneration rate
    k <- exp(dmu_log(x[i,1:d]) - dtarg_log(x[i,1:d]))

    # Time of global move
    if (k > 10^-300){
      t_global <- rexp(1, k)
    } else {
      # (Avoids NAs being produced)
      t_global <- Inf
    }

    if (t_local < t_global){
      # Record holding time
      x[i,d+1] <- t_local
      # Make a local move
      x[i+1,1:d] <- mh_kern(x[i,1:d])
      i <- i+1
    } else {
      # Record holding time
      x[i,d+1] <- t_global
      # Tour ends: regeneration occurs
      regen_yet_to_occur <- 0
    }
  }
  x[1:i,]
}

#' Jump Process Restore Sampler when the jump chain of
#' the underlying jump process has a reversible transition kernel.
#'
#' Assumes the holding rate is 1.
#' Assumes C * (normalising constant of target) = 1.
#'
#' @param ntours Number of tours to simulate
#' @param dtarg_log Function to evaluate the log target density
#' @param dmu_log Function to evaluate the log regeneration distribution
#' @param rmu Function to simulate n times from the regeneration distribution
#' @param mh_kern Metrpolis-Hastings kernel.
#' Input: current state. Output next state.
#' @return A matrix with first d rows (d the dimension of the target density)
#' the jump chain and last row the holding time for the jump process.
jprs <- function(ntours, dtarg_log, dmu_log, rmu, mh_kern){

  jp <- replicate(ntours, jprs_tour(dtarg_log, dmu_log, rmu, mh_kern))
  do.call(rbind, jp)
}
