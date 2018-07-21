#' S4 permutation test class
#' @exportClass prmt
setClass("prmt", representation(permutation_results = "data.frame",
                                original_permutation_results = "list"),
         validity = function(object) {
           errors <- c()
           if (!statistic_names_valid(colnames(object@permutation_results))) {
             error <- "One or more permutation statistic name(s) missing"
             errors <- c(errors, error)
           }
           
           if (!statistic_names_valid(names(object@original_permutation_results)) ||
               ncol(object@permutation_results) != length(object@original_permutation_results) ||
               any(colnames(object@permutation_results) != names(object@original_permutation_results))) {
             error <- "One or more original permutation statistic name(s) missing"
             errors <- c(errors, error)
           }
           
           if (length(errors) == 0) TRUE else errors
         })

#' S4 permutation test summary  class
#' @exportClass summary.prmt
setClass("summary.prmt", representation(summary = "data.frame",
                                        permutations = "integer",
                                        level = "numeric",
                                        alternative = "character"),
         validity = function(object) {
           errors <- c()
           
           required_summary_columns <- c("statistic", "ci_lower", "ci_upper", "actual", "p_value")
           if (is.null(colnames(object@summary)) ||
               length(setdiff(required_summary_columns, colnames(object@summary))) != 0 &&
               length(required_summary_columns) != length(colnames(object@summary))) {
             errors <- c(errors,
                         paste0("Summary column names do not match expectation of ", paste(required_summary_columns, collapse = ',')))
           }
           
           if (any(object@summary$p_value < 0 | object@summary$p_value > 1)) {
             errors <- c(errors,
                         "Summary contains a p_value that does not represent a valid probability")
           }
           
           if (!(object@level <= 1 && object@level >= 0)) {
             errors <- c(errors,
                         "Significance level must fall in [0,1]")
           }
           
           valid_alternatives <- c('less', 'greater', 'two.sided')
           if (!(object@alternative %in% valid_alternatives)) {
             errors <- c(errors,
                         paste0("Alternative must be one of ", paste(valid_alternatives, collapse = ',')))
           }
           
           if (length(errors) == 0) TRUE else errors
         })

#' Perform a full permutation test
#'
#' This function should only be used for academic purposes or small sample sizes (time grows factorially with number of observations)
#' If the "parallel" pacakge is installed, use the global mc.cores option to control level of parallelism
#'
#' @param data a data frame to draw statistics from
#' @param statistics a function or vector of functions accepting a data frame, returning a scalar statistic
#' @param group_var a character, singular column name in the data representing groups.
#' @param statistic_names the names of statistics that are computed. must be specified if the supplied statistics do not have names
#' @param trials the number of Monte Carlo trials to perform
#' 
#' @author kholub
#' @examples
#' library(dplyr)
#' data(mpg)
#' test_data <- mpg %>% 
#'   filter(class %in% c('suv', 'compact')) %>%
#'   sample_n(7) # this just for computational feasibility of illustrative example
#' prm_res <- prm_test(test_data, c(median_difference = function(df) {
#'                                                                      group_medians <- df %>%
#'                                                                        group_by(class) %>%
#'                                                                        summarize(
#'                                                                          median = median(hwy))
#'                                                                      as.numeric(group_medians[group_medians$class == 'suv', 'median']) -
#'                                                                      as.numeric(group_medians[group_medians$class == 'compact', 'median'])
#'                                                                   }), 
#'                     'class')
#'
#' @export prm_test
prm_test <- function(data, statistics, group_var,
                     statistic_names = names(statistics),
                     max_observations = 10) {
  if (nrow(data) > max_observations) {
    stop('Attempting to perform a permutation test of', as.character(max_observations), 'observations with max configured to', as.character(max_observations), 
         '. Consider Monte Carlo randomization (approximate) alternatives.')
  }
  
  if (is.function(statistics)) {
    statistics <- c(statistics)
  }
  
  if (length(statistics) < 1) {
    stop('Must specify one or more statistics to test')
  }
  
  if (!statistic_names_valid(statistic_names) ||
      length(statistics) != length(statistic_names)) {
    stop('One or more statistic name(s) missing')
  }
  
  apply_method <- get_apply_method()
  
  # TODO calculate permutations dynamically so we aren'y memory bound
  all_permutations <- generate_all_permutation_indices(nrow(data))
  
  results <- apply_method(1:nrow(all_permutations), function(trial_ix) {
    permuted_ix <- sample(nrow(data), nrow(data))
    permuted_data <- data
    permuted_data[[group_var]] <- data[[group_var]][permuted_ix]
    
    results_row <- lapply(statistics, function(stat){
      stat(permuted_data)
    })
    
    results_row
  })
  
  # better alternatives, but keeps it base R only
  results <- data.frame(do.call(rbind, results))
  colnames(results) <- statistic_names
  
  original_permutation_results <- lapply(statistics, function(stat) {
    stat(data)
  })
  names(original_permutation_results) <- statistic_names
  
  new('prmt',
      permutation_results = results,
      original_permutation_results = original_permutation_results)
}

#' Visualize observed statistics against null distributions
#'
#' @param object an prmt object
#' @param nrows the number of rows in the resulting graphic (each entry corresponds to a single statistic)
#' @param ncols the number of columns in the resulting graphic
#'
#' @author kholub
#'
#' @examples
#' data(mpg)
#' sample1 <- mpg %>% filter(class == 'suv')
#' sample2 <- mpg %>% filter(class == 'compact')
#' mcr_results <- mcrd_test(sample1, sample2, median = median(hwy),
#'                                           mean = mean(hwy))
#' plot(mcr_results)
#'
#' @export plot.prmt
plot.prmt <- function(x, ..., nrows = 1, ncols = NULL) {
  ncols <- if (is.null(ncols)) ncol(x@permutation_results) else ncols
  
  previous_state <- par(mfrow = c(nrows, ncols))
  tryCatch({
    for (stat_name in colnames(x@permutation_results)) {
      stat_permuted <- unlist(x@permutation_results[[stat_name]])
      stat_observed <- unlist(x@original_permutation_results[[stat_name]])
      
      x_limits = c(min(stat_permuted, stat_observed),
                   max(stat_permuted, stat_observed))
      
      hist(stat_permuted, main = stat_name,
           xlab = "Statistic under permutation",
           xlim = x_limits)
      abline(v = stat_observed, col = 'blue')
    }
  },
  finally = { par(previous_state) })
}

setMethod("plot", signature("prmt"), plot.prmt)

#' Compute empirical confidence intervals for statistics under permutation testing
#'
#' @param object an object of class prmt
#' @param parm a vector of statistic names to produce intervals for
#' @param level the width of confidence intervals
#'
#' @author kholub
#' 
#' @examples
#' data(mpg)
#' sample1 <- mpg %>% filter(class == 'suv')
#' sample2 <- mpg %>% filter(class == 'compact')
#' mcr_results <- mcrd_test(sample1, sample2, median = median(hwy),
#'                                            mean = mean(hwy))
#' confint(mcr_results, parm = c("mean", "median"))
#'
#' @export confint.prmt
confint.prmt <- function(object, parm, level = .95, ...) {
  stopifnot(level >= 0 && level <= 1)
  
  summary_df <- as.data.frame(summary(object, level = level))
  
  cis<- summary_df[summary_df$statistic %in% parm, c('ci_lower', 'ci_upper')]
  percentiles <- as.character(100 * c((1 - level) / 2, level + (1 - level) / 2))
  colnames(cis) <- sapply(percentiles, function(x){ paste0(x,'%') })
  rownames(cis) <- parm
  
  as.matrix(cis)
}

setMethod("confint", signature("prmt"), confint.prmt)

#' Print a permutation test
#' 
#' @author kholub
#'
#' @export print.prmt
print.prmt <- function(x, ...) { 
  cat("Monte Carlo randomization test of", as.character(nrow(x@permutation_results)), 
      " trials on the following statistics:", paste0(colnames(x@permutation_results), collapse = ','))
}

setMethod("print", signature(x = "prmt"), print.prmt)
setMethod("show", signature(object = "prmt"), function(object) { print.prmt(object) })

#' Summarize the results of a permutation test
#'
#' Yields p values & CIs against null (empirical) distribution of statistics.
#' For two sided hypotheses, reported p-value is the minimum one sided p-value with a Bonferroni correction.
#'
#' @param object a prmt object to be summarized
#' @param level the width of empirical confidence intervals
#' @param alternative the alternative used for tests of significance. one of 'less', 'greater', 'two.sided'
#'
#' @author kholub
#'
#' @examples
#' data(mpg)
#' sample1 <- mpg %>% filter(class == 'suv')
#' sample2 <- mpg %>% filter(class == 'compact')
#' mcr_results <- mcrd_test(sample1, sample2, mean = mean(hwy))
#' summary(mcr_results)
#'
#' @export summary.prmt
summary.prmt <- function(object, level = .95, alternative = 'two.sided', ...) {
  stopifnot(level <= 1 && level >= 0)
  
  total_summary <- data.frame()
  for (statistic_name in colnames(object@permutation_results)) {
    permuted_differences <- unlist(object@permutation_results[[statistic_name]])
    observed_difference <- object@original_permutation_results[[statistic_name]]
    cdf <- ecdf(permuted_differences)
    
    lower_tail <- cdf(observed_difference)
    if (alternative == 'less') {
      empirical_p <- lower_tail
      empirical_ci <- c(-Inf, quantile(permuted_differences, level))
    }
    else if (alternative == 'greater') {
      empirical_p <- 1 - lower_tail
      empirical_ci <- c(quantile(permuted_differences, 1 - level), Inf)
    }
    else if (alternative == 'two.sided') {
      # by selecting the minimum, we are implicitly performing 2 tests, so we apply a Bonferroni correction
      empirical_p <- 2 * min(1 - lower_tail, lower_tail)
      empirical_ci <- quantile(permuted_differences, c((1 - level) / 2, level + (1 - level) / 2))
    }
    else {
      stop('Unrecognized "alternative" argument. Must be one of "less", "greater", "two.sided"')
    }
    total_summary <- rbind(total_summary, list(statistic = statistic_name,
                                               ci_lower = empirical_ci[1],
                                               ci_upper = empirical_ci[2],
                                               actual = observed_difference,
                                               p_value = empirical_p),
                           stringsAsFactors = FALSE)
  }
  
  rownames(total_summary) <- NULL
  
  new('summary.prmt', 
      summary = total_summary,
      permutations = nrow(object@permutation_results),
      level = level,
      alternative = alternative)
}

setMethod("summary", signature("prmt"), summary.prmt)

#' Access permutation test summary data
#'
#' @author kholub
#'
#' @export as.data.frame.summary.prmt
as.data.frame.summary.prmt <- function(x, ...) { 
  x@summary
}

setMethod("as.data.frame", signature("summary.prmt"), as.data.frame.summary.prmt)

#' Print a permutation test summary
#' 
#' @author kholub
#'
#' @export print.summary.prmt
print.summary.prmt <- function(x, ...) { 
  cat("Permutation test of", x@permutations, "permutations against alternative of", x@alternative, "at significance", as.character(x@level), '\n')
  print(x@summary, row.names = FALSE)
}

setMethod("print", signature(x = "summary.prmt"), print.summary.prmt)
setMethod("show", signature(object = "summary.prmt"), function(object) { print.summary.prmt(object) })