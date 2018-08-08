#' S4 permutation test power class
#' @exportClass prm_power
setClass("prm_power", representation(experiments = "list",
                                     sample_sizes = "data.frame"), # allow for multiple groups in the case of imbalanced designs
         validity = function(object) {
           errors <- c()
           
           if (length(object@experiments) < 1) {
             error <- 'Must simulate at least one experiment to determine power'
             errors <- c(errors, error)
           }
           
           if (!is(object@experiments[[1]], 'prmt')) {
             error <- 'Experiments must be instances of prmt'
             errors <- c(errors, error)
           }
           
           statistics <- colnames(object@experiments[[1]]@permutation_results)
           if (!statistic_names_valid(statistics)) {
             error <- "One or more statistic name(s) missing from permutation test results"
             errors <- c(errors, error)
           }
           
           if (nrow(object@sample_sizes) < 1) {
             error <- "Must provide valid sample sizes"
             errors <- c(errors, error)
           }
           
           if (length(errors) == 0) TRUE else errors
         })

#' S4 randomization test power class
#' @exportClass mcr_power
setClass("mcr_power", contains = "prm_power")

#' S4 permutation test power summary class
#' @exportClass summary.prm_power
setClass("summary.prm_power", representation(statistic_summary = "data.frame",
                                             sample_sizes = "data.frame", # allow for multiple groups in the case of imbalanced designs
                                             alternative = 'character',
                                             alpha = 'numeric',
                                             experiments = 'integer'), 
         validity = function(object) {
           errors <- c()
           
           required_summary_columns <- c('statistic', 'power', 'average_effect')
           if ( length(setdiff(required_summary_columns, colnames(object@statistic_summary))) != 0 &&
                length(required_summary_columns) != length(colnames(object@statistic_summary))) {
             error <- 'Missing required statistic summary columns'
             errors <- c(errors, error)
           }
           
           if (!statistic_names_valid(object@statistic_summary$statistic)) {
             error <- "One or more statistic name(s) missing from power summary"
             errors <- c(errors, error)
           }
           
           if (nrow(object@sample_sizes) < 1) {
             error <- "Must provide valid sample sizes"
             errors <- c(errors, error)
           }
           
           if (length(errors) == 0) TRUE else errors
         })

#' Perform a power analysis of a randomization test
#'
#' @param generate_data a function returning a dataframe matching structure expected by statistics
#' @param statistics a list of functions to evaluate data produced by generate_data
#' @param group_var a character column name to split the data into groups
#' @param test_trials the number of MC trials to be performed per experiment
#' @param power_trials the number of experiments to be performed in the power calculation
#'
#' @author kholub
#' 
#' @examples 
#' gen_data <- function(){ data.frame(x = c(rnorm(50), rnorm(50, 1)), class = rep(c('a', 'b'), each = 50)) }
#' mcrp <- mcr_power(gen_data, c(mean = function(df){ mean(df[df$class == 'a', 'x']) - mean(df[df$class == 'b', 'x']) }, 
#'                               median = function(df){ median(df[df$class == 'a', 'x']) - median(df[df$class == 'b', 'x']) }),
#'     'class')
#' mcrp
#' 
#' @export mcr_power
mcr_power <- function(generate_data, statistics, group_var,
                      test_trials = 1e2,
                      power_trials = 1e2) {
  stopifnot(power_trials > 0)
  stopifnot(test_trials > 0)
  
  # since mcr_test can run in parallel, always run serially to avoid nested parallelization
  experiments <- lapply(1:power_trials, function(trial_ix) {
    mcr_test(generate_data(), statistics, group_var, 
             trials = test_trials)
  })
  
  group_sizes <- as.data.frame(table(generate_data()[,group_var]))
  names(group_sizes) <- c(group_var, 'size')
   
  new('mcr_power',
      experiments = experiments,
      sample_sizes = group_sizes)
}

#' Perform a power analysis of a randomization test
#' 
#' Note this is a convenience wrapper for \code{\link{mcr_power}} in the case of two sample differences
#' 
#' @param generate_sample1 a function returning a dataframe matching structure expected by statistics
#' @param generate_sample2 ''
#' @param ... expressions to evaluate data produced by sample generating functions
#' @param test_trials the number of MC trials to be performed per experiment
#' @param power_trials the number of experiments to be performed in the power calculation
#'
#' @author kholub
#' @examples 
#' gen_data_s1 <- function(){ data.frame(x = rnorm(50), class = 'a') }
#' gen_data_s2 <- function(){ data.frame(x = rnorm(50, 1), class = 'b') }
#' mcrp <- mcrd_power(gen_data_s1, gen_data_s2, mean = mean(x), median = median(x), 
#'                    test_trials = 1e2)
#' mcrp
#'
#' @export mcrd_power
mcrd_power <- function(generate_sample1, generate_sample2, ...,
                       test_trials = 1e2,
                       power_trials = 1e2) {
  stopifnot(power_trials > 0)
  stopifnot(test_trials > 0)
  
  # since mcr_test can run in parallel, always run serially here to avoid nested parallelization
  experiments <- lapply(1:power_trials, function(trial_ix) {
    mcrd_test(generate_sample1(), generate_sample2(), ..., trials = test_trials)
  })
  
  group_sizes <- data.frame(sample = c('1', '2'), 
                            size = c(nrow(generate_sample1()), nrow(generate_sample2())))
  
  new('mcr_power',
      experiments = experiments,
      sample_sizes = group_sizes)
}

#' Visualize the power analysis of a randomization test
#' 
#' If alpha is left unspecified, the distribution of p-values is visualized. If alpha is specified, power is plotted as a function of alpha
#'
#' @param x a prm_power object to be plotted
#' @param statistic the name of the statistic to be visualized
#' @param alternative the alternative hypotheses for p-value determination. one of "less", "greater", "two.sided"
#' @param alpha the numeric significance level(s) [0-1] used in inference - determines the FPR
#'
#' @author kholub
#' @examples
#' gen_data <- function(){ data.frame(x = c(rnorm(50), rnorm(50, 1)), class = rep(c('a', 'b'), each = 50)) }
#' mcrp <- mcr_power(gen_data, c(mean = function(df){ mean(df[df$class == 'a', 'x']) - mean(df[df$class == 'b', 'x']) }, 
#'                               median = function(df){ median(df[df$class == 'a', 'x']) - median(df[df$class == 'b', 'x']) }),
#'     'class')
#' plot(mcrp, statistic = 'median', alternative = 'two.sided')
#' plot(mcrp, statistic = 'median', alternative = 'two.sided', alpha = NULL)
#'
#' @export plot.prm_power
plot.prm_power <- function(x, ..., statistic, alternative, 
                           alpha = seq(.01, .25, by = .01)) {
  stopifnot(statistic %in% colnames(x@experiments[[1]]@permutation_results))
  stopifnot(is.null(alpha) || is.na(alpha) || length(alpha) > 0)
  
  if (is.null(alpha) || is.na(alpha)) {
    p <- sapply(x@experiments, function(experiment){
      s <- as.data.frame(summary(experiment, alternative = alternative))
      s[s$statistic == statistic, 'p_value']
    })
    
    hist(p, 
         main = paste0('p-values for ', statistic),
         xlab = 'p-value',
         freq = TRUE)
  }
  else {
    power <- sapply(alpha, function(sig) {
      s <- as.data.frame(summary(x, alternative = alternative, alpha = sig))
      s[s$statistic == statistic, 'power']
    })
    
    plot(alpha, power, 
         main = paste0('Power analysis of ', statistic),
         xlab = 'alpha',
         ylab = 'power')
  }
}

setMethod("plot", signature("prm_power"), plot.prm_power)

#' Print a permutation test power analysis
#' 
#' @param x a prm_power object to be printed
#' 
#' @author kholub
#' 
#' @examples
#' gen_data <- function(){ data.frame(x = c(rnorm(50), rnorm(50, 1)), class = rep(c('a', 'b'), each = 50)) }
#' mcrp <- mcr_power(gen_data, c(mean = function(df){ mean(df[df$class == 'a', 'x']) - mean(df[df$class == 'b', 'x']) }, 
#'                               median = function(df){ median(df[df$class == 'a', 'x']) - median(df[df$class == 'b', 'x']) }),
#'     'class')
#' print(mcrp)   
#'
#' @export print.prm_power
print.prm_power <- function(x, ...) { 
  cat("Power analysis of a permutation test of", as.character(length(x@experiments)), 
      "experiments on the following statistics:", paste0(colnames(x@experiments[[1]]@permutation_results), collapse = ','))
}

setMethod("print", signature(x = "prm_power"), print.prm_power)
setMethod("show", signature(object = "prm_power"), function(object) { print.prm_power(object) })

#' Summarize a power analysis of a randomization test
#'
#' @param object a prm_power object to be summarized
#' @param alternative the alternative hypotheses for p-value determination. one of "less", "greater", "two.sided"
#' @param alpha the numeric significance level(s) [0-1] used in inference - determines the FPR
#'
#' @author kholub
#' 
#' @examples 
#' gen_data <- function(){ data.frame(x = c(rnorm(50), rnorm(50, 1)), class = rep(c('a', 'b'), each = 50)) }
#' mcrp <- mcr_power(gen_data, c(mean = function(df){ mean(df[df$class == 'a', 'x']) - mean(df[df$class == 'b', 'x']) }, 
#'                               median = function(df){ median(df[df$class == 'a', 'x']) - median(df[df$class == 'b', 'x']) }),
#'     'class')
#' summary(mcrp, .05, 'two.sided')
#'
#'@export summary.prm_power
summary.prm_power <- function(object, alpha, alternative, ...) {
  simulation_list <- lapply(object@experiments, function(experiment) {
    prm_summary <- as.data.frame(summary(experiment, alternative = alternative))
    
    as.data.frame(prm_summary[, c('p_value', 'actual', 'statistic')])
  })
  
  simulation_results <- as.data.frame(do.call(rbind, simulation_list))
  
  summary_list <- lapply(unique(simulation_results$statistic), function(statistic) {
    statistic_simulation <- simulation_results[simulation_results$statistic == statistic,]
    
    list(statistic = statistic,
         power = sum(statistic_simulation$p_value <= alpha) / nrow(statistic_simulation),
         average_effect = mean(statistic_simulation$actual))
  })
  
  statistic_summary <- as.data.frame(do.call(rbind, summary_list))
  
  new('summary.prm_power',
      statistic_summary = statistic_summary,
      sample_sizes = object@sample_sizes,
      alternative = alternative,
      alpha = alpha,
      experiments = length(object@experiments))
}

setMethod("summary", signature("prm_power"), summary.prm_power)

#' Access a permutation test power analysis summarization
#'
#' @param x a prm_power object 
#'
#' @author kholub
#' 
#' @examples
#' gen_data <- function(){ data.frame(x = c(rnorm(50), rnorm(50, 1)), class = rep(c('a', 'b'), each = 50)) }
#' mcrp <- mcr_power(gen_data, c(mean = function(df){ mean(df[df$class == 'a', 'x']) - mean(df[df$class == 'b', 'x']) }, 
#'                               median = function(df){ median(df[df$class == 'a', 'x']) - median(df[df$class == 'b', 'x']) }),
#'     'class')
#' s <- summary(mcrp, .05, 'two.sided')
#' as.data.frame(s)
#'
#' @export as.data.frame.summary.prm_power
as.data.frame.summary.prm_power <- function(x, ...) {
  x@statistic_summary
}

setMethod("as.data.frame", signature("summary.prm_power"), as.data.frame.summary.prm_power)

#' Print a power analysis summary of a permutation test
#'
#' @param x a summary.prm_power object to be printed
#'
#' @author kholub
#' 
#' @examples
#' gen_data <- function(){ data.frame(x = c(rnorm(50), rnorm(50, 1)), class = rep(c('a', 'b'), each = 50)) }
#' mcrp <- mcr_power(gen_data, c(mean = function(df){ mean(df[df$class == 'a', 'x']) - mean(df[df$class == 'b', 'x']) }, 
#'                               median = function(df){ median(df[df$class == 'a', 'x']) - median(df[df$class == 'b', 'x']) }),
#'     'class')
#' s <- summary(mcrp, .05, 'two.sided')
#' print(s)
#'
#'@export print.summary.prm_power
print.summary.prm_power <- function(x, ...) {
  cat("Power analysis of", x@experiments, "experiments with alternative of", x@alternative, "at significance", as.character(x@alpha), '\n')
  cat("Group sizes:\n")
  print(x@sample_sizes, row.names = FALSE)
  print(x@statistic_summary, row.names = FALSE)
}

setMethod("print", signature(x = "summary.prm_power"), print.summary.prm_power)
setMethod("show", signature(object = "summary.prm_power"), function(object) { print.summary.prm_power(object) })
