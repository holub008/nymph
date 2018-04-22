#' S4 Monte Carlo randomization test class
#' @exportClass mcrt
setClass("mcrt", contains = "prmt")
# preferring inheritance to delegation since for boilerplate reduction

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
  
  new('prmt',
      permutation_results = results,
      original_permutation_results = original_permutation_results)
}
