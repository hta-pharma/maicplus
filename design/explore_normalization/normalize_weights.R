normalize_weights <- function(weights, method) {
  n <- length(weights)
  
  # Handle edge case of zero weights
  if (sum(weights) == 0) {
    return(rep(1, n))
  }
  
  switch(method,
         # Original weights (no normalization)
         "OW" = weights,
         
         # Sum to one (standard normalization)
         # Manuscript: w̃_i = w_i / Σw_j
         "SW1" = weights / sum(weights),
         
         # Sum to N (preserve sample size interpretation)
         # Manuscript: w̃_i = n × w_i / Σw_j
         "SWN" = weights * n / sum(weights),
         
         # Sum to ESS (balance between SW1 and SWN)
         # Manuscript: w̃_i = ESS(w) × w_i / Σw_j
         "SWESS" = {
           ess <- sum(weights)^2 / sum(weights^2)
           weights * ess / sum(weights)
         },
         
         stop("Unknown normalization method")
  )
}