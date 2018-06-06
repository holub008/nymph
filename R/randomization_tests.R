#' S4 Monte Carlo randomization test class
#' @exportClass mcrt
setClass("mcrt", contains = "prmt")
# preferring inheritance to delegation since for boilerplate reduction

#' Perform a Monte Carlo randomization test
#' 
#' Simulate a permutation test of arbitrary statistics using a simple randomization (i.e. non-clustered) scheme
#' Note that \code{\link{mcrd_test}} and \code{\link{mcr_anova}} convenience functions cover common use cases and should be preferred for interactive, non-programmatic use
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
#' test_data <- mpg %>% filter(class %in% c('suv', 'compact'))
#' mcr_res <- mcr_test(test_data, c(median_difference = function(df) {
#'                                                                       group_medians <- df %>%
#'                                                                         group_by(class) %>%
#'                                                                         summarize(
#'                                                                           median = median(hwy))
#'                                                                       as.numeric(group_medians[group_medians$class == 'suv', 'median']) -
#'                                                                          as.numeric(group_medians[group_medians$class == 'compact', 'median'])
#'                                                                      }),
#'                       'class')
#'
#' @export mcr_test
mcr_test <- function(data, statistics, group_var,
                     statistic_names = names(statistics),
                     trials = 1e3) {
  if (trials < 1) {
    stop('Must specify a positive number of trials to run')
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
  
  results <- apply_method(1:trials, function(trial_ix) {
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
  
  new('mcrt',
      permutation_results = results,
      original_permutation_results = original_permutation_results)
}

#' Perform a two sample Monte Carlo randomization test for difference in arbitrary statistics
#'
#' Randomization testing for assessing significance of difference in arbitrary test statistics across two samples.
#' If the "parallel" pacakge is installed, use the global mc.cores option to control level of parallelism
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
#'                                         median = median(hwy))
#'
#' @export mcrd_test
mcrd_test <- function(sample1, 
                      sample2,
                      ...,
                      trials = 1e3) {
  
  # TODO delegate the meat of this to mcr_test
  
  if (trials < 1) {
    stop('Must specify a positive number of trials to run')
  }
  
  # note that the call object includes a sentinel first value - we ignore it
  statistics_call <- substitute(list(...))[-1]
  statistic_names <- names(statistics_call)
  eval_env <- parent.frame()
  
  if (length(statistic_names) < 1) {
    stop('Must specify one or more statistic(s) to test')
  }
  else if (is.null(statistic_names) || any(statistic_names == '')) {
    stop('One or more statistic name(s) missing')
  }
  
  total_sample <- rbind(sample1, sample2)
  
  apply_method <- get_apply_method()
  
  results <- apply_method(1:trials, function(trial_ix) {
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
