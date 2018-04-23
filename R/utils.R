# location for shared functionality

get_apply_method <- function() {
  if (requireNamespace("parallel", quietly = TRUE)) {
    return(parallel::mclapply)
  }
  else {
    return(lapply)
  }
}

# combinat::permn accomplishes something similar, but is a needless dependency due to simplicity
# note that this would ideally be inlined to avoid storing results in memory, but is a premature optimization at current
# alternatively, a c implementation may be appropriate
generate_all_permutation_indices <- function(n) {
  stopifnot(n > 0)
  
  if (n == 1) {
    return(matrix(1, nrow = 1))
  }
  else {
    sub_problem_perms <- generate_all_permutation_indices(n - 1)
    current_problem_perms <- matrix(nrow = nrow(sub_problem_perms) * n, ncol = n)
    # this clever (dare I say obfuscated?) algo stolen from SO
    for (candidate_ix in 1:n) {
      current_problem_perms[(candidate_ix - 1) * nrow(sub_problem_perms) + (1:nrow(sub_problem_perms)), ] <- 
        cbind(candidate_ix, sub_problem_perms + (sub_problem_perms >= candidate_ix))
    }
    
    return(current_problem_perms)
  }
}

statistic_names_valid <- function(statistic_names) {
  # sneakily relies on short circuit evaluation to avoid a warning on NULL input
  length(statistic_names) > 0 && 
    !any(is.na(statistic_names)) &&
    !any(statistic_names == '')
}
