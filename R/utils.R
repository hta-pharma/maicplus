# Restore the RNG back to a previous state using the global .Random.seed
set_random_seed <- function(old_seed) {
  if (is.null(old_seed)) {
    rm(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    assign(".Random.seed", value = old_seed, envir = globalenv(), inherits = FALSE)
  }
}
