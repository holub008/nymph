#' S4 bootstrap class
#' @exportClass mcrt
setClass("mcrt", representation(permutation_results = "data.frame",
                                original_permutation_results = "list"),
         validity = function(object) {
           errors <- c()
           if (any(is.na(colnames(object@permutation_results))) || # sneakily relies on short circuit evaluation
               is.null(colnames(object@permutation_results)) ||
               any(colnames(object@permutation_results) == '')) {
             error <- "One or more bootstrapped statistic name(s) missing"
             errors <- c(errors, error)
           }
           
           if (any(is.na(names(object@original_permutation_results))) ||
               is.null(names(object@original_permutation_results)) ||
               any(colnames(object@permutation_results) != names(object@original_permutation_results))) {
             error <- "One or more bootstrapped original permutation statistic name(s) missing"
             errors <- c(errors, error)
           }
           
           if (length(errors) == 0) TRUE else errors
         })

#' S4 bootstrap class
#' @exportClass summary.mcrt
setClass("summary.mcrt", representation(summary = "data.frame",
                                trials = "integer",
                                significance_level = "numeric",
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
           
           if (!(significance_level <= 1 && significance_level >= 0)) {
             errors <- c(errors,
                         "Significance level must fall in [0,1]")
           }
           
           valid_alternatives <- c('less', 'greater', 'two.sided')
           if (!(alternative %in% valid_alternatives)) {
             errors <- c(errors,
                         paste0("Alternative must be one of ", paste(valid_alternatives, collapse = ',')))
           }
           
           if (length(errors) == 0) TRUE else errors
         })

#' Perform a two sample Monte Carlo randomization test for difference in arbitrary statistics
#'
#' Randomization testing for assessing significance of difference in arbitrary test statistics across two samples
#' Monte Carlo sampling is used to make the test computationally feasible
#'
#' @param sample1 a data.frame representing observations from a single sample or treatment
#' @param sample2 ''
#' @param ... named statistics to be computed on the 
#' @param trials the number of Monte Carlo trials to perform
#' 
#' @return a Monte Carlo Randomization Test object containing the results of the test
#'
#' @author kholub
#' @examples
#' data(mpg)
#' sample1 <- mpg %>% filter(class == 'suv')
#' sample2 <- mpg %>% filter(class == 'compact')
#' mcrd_res <- mcrd_test(sample1, sample2, mean = mean(hwy),
#'                                            median)
#'
#' @importFrom parallel mclapply
#'
#' @export mcrd_test
mcrd_test <- function(sample1, 
                      sample2,
                      ...,
                      trials = 1e3) {
  # note that the call object includes a sentinel first value - we ignore it
  statistics_call <- substitute(list(...))[-1]
  statistic_names <- names(statistics_call)
  eval_env <- parent.frame()
  
  if (length(statistic_names) < 1) {
    stop('Must specify one or more statistics to bootstrap')
  }
  else if (is.null(statistic_names) || any(statistic_names == '')) {
    stop('One or more statistic(s) specified without a name - all statistics must have a name')
  }
  
  total_sample <- rbind(sample1, sample2)
  
  results <- mclapply(1:trials, function(trial_ix) {
    resample1_ix <- sample(nrow(total_sample), nrow(sample1))
    resample1 <- total_sample[resample1_ix,]
    resample2 <- total_sample[-resample1_ix,]
    
      results_row <- lapply(unname(statistics_call), function(stat){
        sample1_stat <- eval(stat, resample1, eval_env)
        sample2_stat <- eval(stat, resample2, eval_env)
        
        sample1_stat - sample2_stat
      })
  })
  
  # better alternatives, but keeps it base R only
  results <- data.frame(do.call(rbind, results))
  
  colnames(results) <- statistic_names
  
  original_permutation_results <- lapply(unname(statistics_call), function(stat) {
    sample1_stat <- eval(stat, sample1, eval_env)
    sample2_stat <- eval(stat, sample2, eval_env)
      
    sample1_stat - sample2_stat
  })
  names(original_permutation_results) <- statistic_names
  
  new('mcrt',
      permutation_results = results,
      original_permutation_results = original_permutation_results)
}

#' Visualize the distribution of statistic comparisons between two samples
#'
#' @param object an mcrt object
#' @param nrows the number of rows in the resulting graphic (each entry corresponds to a single statistic)
#' @param ncols the number of columns in the resulting graphic
#'
#' @author kholub
#'
#' @examples
#' data(mpg)
#' sample1 <- mpg %>% filter(class == 'suv')
#' sample2 <- mpg %>% filter(class == 'compact')
#' mcr_results <- mcr_test(sample1, sample2, median = median(hwy),
#'                                           mean = mean(hwy))
#' plot(mcr_results)
#'
#' @export plot.mcrt
plot.mcrt <- function(x, ..., nrows = 1, ncols = NULL) {
  ncols <- if (is.null(ncols)) ncol(object@permutation_results) else ncols
  
  previous_state <- par(mfrow = c(nrows, ncols))
  tryCatch({
    for (stat_name in colnames(object@permutation_results)) {
      stat_difference_permuted <- unlist(object@permutation_results[[stat_name]])
      stat_difference_observed <- unlist(object@original_permutation_results[[stat_name]])
      
      x_limits = c(min(stat_difference_permuted, stat_difference_observed),
                   max(stat_difference_permuted, stat_difference_observed))
      
      hist(stat_difference_permuted, main = stat_name,
           xlab = "Statistic differences under randomization",
           xlim = x_limits)
      abline(v = stat_difference_observed, col = 'blue')
    }
  },
  finally = { par(previous_state) })
}

setMethod("plot", signature("mcrt"), plot.mcrt)

#' Produce confidence intervals for statistics under Monte Carlo randomization testing
#'
#' @param object an object of class mcrt
#' @param parm a vector of statistic names to produce intervals for
#' @param level the width of confidence intervals
#'
#' @author kholub
#' 
#' @examples
#' data(mpg)
#' sample1 <- mpg %>% filter(class == 'suv')
#' sample2 <- mpg %>% filter(class == 'compact')
#' mcr_results <- mcr_test(~ hwy, sample1, sample2, median = median(hwy))
#' confint(mcr_results, parm = "median")
#'
#' @export confint.mcrt
confint.mcrt <- function(object, parm, level = .95, ...) {
  summary_df <- as.data.frame(summary(object, empirical_interval_level = level))
  
  cis<- summary_df[summary_df$statistic_name %in% parm, c('ci_lower_bound', 'ci_lower_bound')]
  colnames(cis) <- as.character(c((1 - level) / 2, level + (1 - level) / 2))
  
  as.matrix(cis)
}

setMethod("confint", signature("mcrt"), confint.mcrt)

#' Print a Monte Carlo randomization test
#' 
#' @author kholub
#'
#' @export print.summary.mcrt
print.mcrt <- function(x, ...) { 
  cat("Monte Carlo randomization test of", as.character(x@trials), " trials on the following statistics:", paste0(x@permutation_results, collapse = ','))
  cat(x@summary)
}

setMethod("print", signature(x = "mcrt"), print.mcrt)
setMethod("show", signature(object = "mcrt"), function(object) { print.mcrt(object) })

#' Summarize the results of a Monte Carlo randomization test
#'
#' Yields p values & CIs against null (empirical) distribution of statistic differences.
#' For two sided hypotheses, reported p-value is the minimum one sided p-value with a Bonferroni correction.
#'
#' @param object an mcrt object to be summarized
#' @param significance_level the significance level for empirical confidence intervals of statistic differences
#' @param alternative the alternative used for tests of significance of statistic differences. one of 'less', 'greater', 'two.sided'
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
#' @importFrom dplyr bind_rows
#'
#' @export summary.mcrt
summary.mcrt <- function(object, significance_level = .05, alternative = 'two.sided', ...) {
  stopifnot(significance_level <= 1 && significance_level >= 0)
  
  total_summary <- data.frame()
  for (statistic_name in colnames(object@permutation_results)) {
    permuted_differences <- unlist(object@permutation_results[[statistic_name]])
    observed_difference <- object@original_permutation_results[[statistic_name]]
    cdf <- ecdf(permuted_differences)
    
    lower_tail <- cdf(observed_difference)
    if (alternative == 'less') {
      empirical_p <- lower_tail
      empirical_ci <- c(-Inf, quantile(permuted_differences, significance_level))
    }
    else if (alternative == 'greater') {
      empirical_p <- 1 - lower_tail
      empirical_ci <- c(quantile(permuted_differences, 1 - significance_level), Inf)
    }
    else if (alternative == 'two.sided') {
      # by selecting the minimum, we are implicitly performing 2 tests, so we apply a Bonferroni correction
      empirical_p <- 2 * min(1 - lower_tail, lower_tail)
      empirical_ci <- quantile(permuted_differences, c(significance_level / 2, 1 - significance_level / 2))
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
  
  new('summary.mcrt', 
      summary = total_summary,
      trials = nrow(object@permutation_results),
      significance_level = significance_level,
      alternative = alternative)
}

setMethod("summary", signature("mcrt"), summary.mcrt)

#' Access Monte Carlo randomization test summary data
#'
#' @author kholub
#'
#' @export as.data.frame.summary.mcrt
as.data.frame.summary.mcrt <- function(x, ...) { 
  x@summary
}

setMethod("as.data.frame", signature("summary.mcrt"), as.data.frame.summary.mcrt)

#' Print a Monte Carlo randomization test summary
#' 
#' @author kholub
#'
#' @export print.summary.mcrt
print.summary.mcrt <- function(x, ...) { 
  cat("Randomization test on", as.character(x@trials), " against alternative of", x@alternative, "at significance", as.character(x@significance_level), '\n')
  cat(x@summary)
}

setMethod("print", signature(x = "summary.mcrt"), print.summary.mcrt)
setMethod("show", signature(object = "summary.mcrt"), function(object) { print.summary.mcrt(object) })
