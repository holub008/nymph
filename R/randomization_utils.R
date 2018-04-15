# location for shared functionality across all tests

get_apply_method <- function() {
  if (requireNamespace("parallel", quietly = TRUE)) {
    return(parallel::mclapply)
  }
  else {
    return(lapply)
  }
}